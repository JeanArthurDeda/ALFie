import bpy
import math
from restir_core import cross, normalize, add, sub, length, dot, mul, mul_s

def generate_cached_mesh(instance):
    def get_material(mat):
        if not mat.use_nodes: return (1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0)
        bsdf = mat.node_tree.nodes.get("Principled BSDF")            
        if not bsdf: return (1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0)
        c = tuple(bsdf.inputs["Base Color"].default_value)
        ems = bsdf.inputs["Emission Color"].default_value[:3]
        s = bsdf.inputs["Emission Strength"].default_value
        ems = (ems[0] * s, ems[1] * s, ems[2] * s)
        return c, ems

    object = instance.object
        
    mesh = object.data
    mats = [get_material(mat.material) for mat in object.material_slots]
    if len(mats) == 0:
        mats = [((1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.0))]
    verts = [tuple(instance.matrix_world @ v.co) for v in mesh.vertices]
    face_inds = [[v for v in p.vertices] for p in mesh.polygons]
    face_nors = [normalize(cross(sub(verts[inds[1]],verts[inds[0]]), sub(verts[inds[2]],verts[inds[0]]))) for inds in face_inds ]
    face_mats = [p.material_index for p in mesh.polygons]

    return (verts, mats, face_inds, face_nors, face_mats)

def generate_cached_area_lights(meshes):
    area_lights = []
    for mesh_index, mesh in enumerate(meshes):
        verts, mats, face_inds, face_nors, face_mats = mesh

        # An mesh can contain multiple area_lights as an area light is defined as a ~planar surface within a mesh. E.g. a cube has 8 lights
        # An area_light contains multiple faces which are joined based on their vicinity and face normal
        # area_light
        #   mesh_index
        #   area
        #   ems_weight - computed from total faces areas and emissive level (which can be set per each face via material)
        #   faces - list of light faces
        # area_light face
        #   area
        #   ems_weight - computed from total triangles area and emissive level which is constant per face
        #   face_index - indices of face within the mesh_index mesh
        #   areas - areas of each triangle inside the face

        # compute the light faces within the mesh
        light_faces = []
        for face_index, mat_index in enumerate(face_mats):
            mat = mats[mat_index]
            c, ems = mat
            ems_weight = length(ems)
            if ems_weight == 0.0: continue
            # compute the area of each triangls
            inds = face_inds[face_index]
            num_tri = len(inds) - 2
            areas = []
            face_area = 0.0
            for i in range(num_tri):
                a = verts[inds[0]]
                b = verts[inds[i + 1]]
                c = verts[inds[i + 2]]

                area = length(cross(sub(b,a), sub(c,a)))
                face_area = face_area + area
                areas.append(area)
            ems_weight = ems_weight * face_area
            light_faces.append((face_area, ems_weight, face_index, areas))
            
        if 0 == len(light_faces): continue
    
        # group adjiacent and coplanar light_faces into lights
        lights = []
        for face in light_faces:
            (face_area, face_ems_weight, face_index, face_areas) = face
            face_nor = face_nors[face_index]
            inds = face_inds[face_index]

            # try to join this face in an existing area light
            joined = False
            for light in lights:
                for light_face in light:
                    (light_face_area, light_face_ems_weight, light_face_index, light_face_areas) = light_face
                    light_face_nor = face_nors[light_face_index]
                    light_inds = face_inds[light_face_index]

                    # normal rejection
                    if dot(face_nor, light_face_nor) < 0.9: continue

                    # check for adjency
                    for i in inds:
                        for light_i in light_inds:
                            if i != light_i: continue
                            light.append(face)
                            joined = True
                            break
                        if joined: break
                    if joined: break
                if joined: break

            if not joined: # add a new light
                lights.append([face])

        # update the area_lights data fron lights
        for faces in lights:
            area = 0.0
            ems_weight = 0.0
            for face in faces:
                (face_area, face_ems_weight, face_index, face_areas) = face
                area = area + face_area
                ems_weight = ems_weight + face_ems_weight
            area_lights.append((mesh_index, area, ems_weight, faces))

    # Last area_lights is the accumulated ems_weights of all lighs
    area_lights_ems_weight = 0.0
    for light in area_lights:
        mesh_index, area, ems_weight, faces = light
        area_lights_ems_weight = area_lights_ems_weight + ems_weight
    area_lights.append(area_lights_ems_weight)

    for l in area_lights: print (l)

    return area_lights
