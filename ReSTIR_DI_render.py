import sys
import os
import bpy
import array
import math
import random
from typing import Tuple
from abc import ABC, abstractmethod
from mathutils import Vector, Matrix
from time import perf_counter
from typing import List

# Add current .blend directory to Python path
blend_dir = os.path.dirname(bpy.data.filepath)
if blend_dir not in sys.path:
    sys.path.append(blend_dir)


# =======
# Helpers
# =======

def get_cosmetic_duration(t):
    ms = int(t * 1000) % 1000
    t = int(t)

    d, t = divmod(t, 86400)
    h, t = divmod(t, 3600)
    m, s = divmod(t, 60)

    parts = [f"{v}{u}" for v, u in [(d,"d"), (h,"h"), (m,"m"), (s,"s")] if v]
    ret = " ".join(parts)
    if len(ret) == 0: return "now"
    return ret

def get_orthonormal_axis(z_axis):
    # generate perp vectors
    x_axis = Vector((1.0, 0.0, 0.0))
    y_axis = z_axis.cross(x_axis)
    if y_axis.length < 0.01:
        x_axis = Vector((0.0, 1.0, 0.0))
        y_axis = z_axis.cross(x_axis)
        
    x_axis = y_axis.cross(z_axis)

    return x_axis, y_axis

def get_weighted_random_index(items, total_weight, get_item_weeight):
    weight = random.uniform(0, total_weight)
    a = 0.0
    for index, i in enumerate(items):
        a += get_item_weeight(i)
        if a >= weight: return index
    return len(items) -1

def interp_barycentric(p0, p1, p2, u, v):
    w = 1.0 - u - v
    return p0*u + p1*v + p2*w


# Area light =================

class AreaLight:
    obj = None
    matrix_world = None
    inv_matrix_world = None
    li = None
    verts = None
    ggx_mat_cache = None
    polygons = []
    area = None
    weight = None

    def __init__(self, instance, ggx_mat_cache):
        self.obj = instance.object
        self.ggx_mat_cache = ggx_mat_cache
        obj = instance.object
        mesh = obj.data
        self.matrix_world = instance.matrix_world.copy()
        self.inv_matrix_world = instance.matrix_world.inverted().copy()
        self.verts = [instance.matrix_world @ v.co for v in mesh.vertices]
        mat = self.ggx_mat_cache.get(obj.material_slots[0].material)
        k_d, k_s, k_r, k_e, k_es = mat
        self.li = k_e * k_es

        self.area = 0.0
        self.polygons = []
        for p in mesh.polygons:
            weight = 0
            tri_weights = []
            for i in range(len(p.vertices) - 2):
                a = self.verts[p.vertices[0]]
                b = self.verts[p.vertices[i + 1]]
                c = self.verts[p.vertices[i + 2]]

                area = (b-a).cross(c-a).length
                w = area * k_es
                tri_weights.append(w)
                weight += w
                self.area += area
            self.polygons.append((weight, tri_weights))
        self.weight = self.area * k_es

        # print (f"Area light {self.obj.name}")
        # print (f"\tArea {self.area}")
        # print (f"\tLi {k_es}")
        # print (f"\tWeight {self.weight}")
        # print (f"\tNum polys {len(self.polygons)}")
        # for p in self.polygons:
        #     weight, tri_weights = p
        #     print (f"\t\tWeight {weight}")
        #     print (f"\t\tNum Tris {len(tri_weights)}")
        #     for w in tri_weights:
        #         print (f"\t\t\ttri weight {w}")

    # returns p, n, mat
    def get_weight_random_point(self):
        poly_index = get_weighted_random_index(self.polygons, self.weight, lambda p: p[0])
        poly_weight, tri_weights = self.polygons[poly_index]
        tri_index = get_weighted_random_index(tri_weights, poly_weight, lambda w : w)
        mesh = self.obj.data
        p = mesh.polygons[poly_index]
        a = self.verts[p.vertices[0]]
        b = self.verts[p.vertices[tri_index + 1]]
        c = self.verts[p.vertices[tri_index + 2]]
        u = random.random()
        v = random.random()
        if u + v > 1.0: 
            u = 1.0 - u
            v = 1.0 - v
        p = interp_barycentric(a, b, c, u, v)
        n = ((b-a).cross(c-a)).normalized()
        return p, n
    
    # returns pdf, wi
    def pdf(self, p, n, p_light, n_light):
        wi = (p_light - p)
        r = wi.length
        wi /= r
        pdf_area_uniform = 1.0 / self.area
        cos_theta_light = math.fabs(n_light.dot(-wi)) # double sided pdf
        if cos_theta_light < 1e-6:
            return 0, wi
        return pdf_area_uniform * r * r / cos_theta_light, wi
    
    # returns hit, pos, nor
    def ray_trace(self, p, wi, d):
        lp = self.inv_matrix_world @ p
        le = self.inv_matrix_world @ (p + wi * d)
        ld = (le-lp).length
        lwi = (le-lp).normalized()
        result, lpos, lnor, face_index = self.obj.ray_cast(lp, lwi, distance=ld)
        if result:
            pos = self.matrix_world @ lpos
            nor = (self.matrix_world.to_3x3() @ lnor).normalized()
            
            return True, pos, nor
        return False, None, None
    

# Samplers ===================

def uniform_sample(n):
        u = random.uniform(0.0, 2.0 * math.pi)
        cos_v = random.uniform(0.0, 1.0)
        v = math.acos(cos_v)

        x = math.sin(v) * math.cos(u)
        y = math.sin(v) * math.sin(u) 
        z = math.cos(v)

        pdf = 1.0 / (2.0 * math.pi)
        
        x_axis, y_axis = get_orthonormal_axis(n)
        wi = x*x_axis + y*y_axis + z*n

        cos_theta = z

        return wi, pdf, cos_theta

def uniform_pdf(p, n, wi):
    return 1.0 / (2.0 * math.pi)

def cosine_sample(n):
        u = random.uniform(0.0, 1.0)
        v = random.uniform(0.0, 2.0 * math.pi)

        r = math.sqrt(u)
        x = r*math.cos(v)
        y = r*math.sin(v)
        z = math.sqrt(1.0-u)

        pdf = z / math.pi

        x_axis, y_axis = get_orthonormal_axis(n)
        wi = x*x_axis + y*y_axis + z*n

        cos_theta = z

        return wi, pdf, cos_theta

def cosine_pdf(p, n, wi):
    cos_theta = wi.dot(n)
    if cos_theta <= 0.0: return 0
    return cos_theta / math.pi

SAMPLER_DISTANCE = 20.0
SAMPLER_BIAS = 0.001
class AreaLightsImportanceSampler:
    area_lights : List[AreaLight] = []
    weight = 0.0
    li = None

    def __init__(self, area_lights : List[AreaLight]):
        self.area_lights = area_lights
        for l in self.area_lights:
            self.weight += l.weight

    def samples(self, p, n, num):
        S = []
        for _ in range(num):
            l = self.area_lights[get_weighted_random_index(self.area_lights, self.weight, lambda l: l.weight)]
            pos, nor = l.get_weight_random_point()
            pdf, wi = l.pdf(p, n, pos, nor)
            pdf *= l.weight / self.weight # adjust pdf based on light selection pdf
            cos_theta = max(0, wi.dot(n))
            S.append((wi, pos, nor, pdf, cos_theta, l))
        return S
    
    def pdf(self, p, n, wi):
        pdf = 0.0
        for l in self.area_lights:
            hit, pos, nor = l.ray_trace(p, wi, SAMPLER_DISTANCE)
            if not hit: continue
            l_pdf, _ = l.pdf(p, n, pos, nor)
            pdf += l_pdf * l.weight / self.weight
        return pdf
    
    def get_name(self):
        return f"LightsImportance({len(self.area_lights)} lights)"
    
# BDRF =========================

def lambert_BDRF(k_d):
    return k_d / math.pi

def phong_BRDF(wi, wo, N):
    k_s = 0.75
    shininess = 32.0
    
    R = 2.0 * N.dot(wi) * N - wi
    
    specular = k_s * math.pow(max(0, R.dot(wo)), shininess) #(shininess + 2)/(2*math.pi)
    
    return specular

def ggx_BDRF(wi: Vector, wo: Vector, N: Vector, F0: Vector, roughness: float):
    """
    Evaluate GGX BRDF for given directions in world space.
    
    Args:
        wi: Light direction (pointing toward light source) - unit vector
        wo: Camera direction (pointing toward camera) - unit vector
        N: Surface normal - unit vector
        F0: Specular reflectance at normal incidence (RGB vector, e.g., (0.04, 0.04, 0.04) for dielectrics)
        roughness: Surface roughness (0 = smooth, 1 = rough)
    
    Returns:
        BRDF value (scalar, same for all RGB channels, multiply by F0 for final color)
    """
    # Ensure all vectors are normalized
    wi = wi.normalized()
    wo = wo.normalized()
    N = N.normalized()
    
    # Half-vector
    H = (wi + wo).normalized()
    
    # Dot products
    NdotWi = max(0.0, wi.dot(N))
    NdotWo = max(0.0, wo.dot(N))
    NdotH = max(0.0, H.dot(N))
    WodotH = max(0.0, wo.dot(H))
    
    # If light or view is below surface, BRDF is 0
    if NdotWi <= 0.0 or NdotWo <= 0.0:
        return 0.0
    
    # GGX distribution term D
    alpha = roughness * roughness  # Roughness squared
    alpha_sq = alpha * alpha
    
    NdotH_sq = NdotH * NdotH
    denom_d = NdotH_sq * (alpha_sq - 1.0) + 1.0
    D = alpha_sq / (math.pi * denom_d * denom_d)
    
    # Smith geometry term G (height-correlated)
    def g1(NdotV):
        """Smith G1 function for GGX."""
        if NdotV <= 0.0:
            return 0.0
        tan_theta_sq = (1.0 - NdotV * NdotV) / (NdotV * NdotV)
        return 2.0 / (1.0 + math.sqrt(1.0 + alpha_sq * tan_theta_sq))
    
    G = g1(NdotWi) * g1(NdotWo)
    
    # Fresnel term F (Schlick approximation)
    F = F0 + (Vector((1.0, 1.0, 1.0)) - F0) * (1.0 - WodotH) ** 5
    
    # GGX BRDF (returns RGB vector)
    denominator = 4.0 * NdotWi * NdotWo
    if denominator == 0.0:
        return Vector((0.0, 0.0, 0.0))
    
    brdf_rgb = (D * G * F) / denominator
    
    # Return as scalar (luminance) or RGB depending on needs
    # For scalar return, use luminance:
    # return 0.2126 * brdf_rgb.x + 0.7152 * brdf_rgb.y + 0.0722 * brdf_rgb.z
    
    # Or return RGB vector
    return brdf_rgb

class GGXMaterialCache:
    mats = None

    def __init__(self):
        self.mats = {}

    def get(self, mat):
        name = mat.name
        if name in self.mats: return self.mats[name]
        bsdf = mat.node_tree.nodes.get("Principled BSDF")            
        k_d = Vector(bsdf.inputs["Base Color"].default_value[:3])
        k_r = bsdf.inputs["Roughness"].default_value
        k_s = Vector(bsdf.inputs["Specular Tint"].default_value[:3])
        k_e = Vector(bsdf.inputs["Emission Color"].default_value[:3])
        k_es = bsdf.inputs["Emission Strength"].default_value
        ggx_mat = (k_d, k_s, k_r, k_e, k_es)
        self.mats[name] = ggx_mat
        # print (f"Caching material {name} as:")
        # print (f"\tk_d {k_d}")
        # print (f"\tk_r {k_r}")
        # print (f"\tk_s {k_s}")
        # print (f"\tk_e {k_e}")
        # print (f"\tk_es {k_es}")
        return ggx_mat
    
# Reservoir
def reservoir_new():
    l = None
    l_pos = None
    l_nor = None
    w_sum = 0
    M = 0
    pdf = 0
    return (l, l_pos, l_nor, w_sum, M, pdf)

def reservoir_update(r, l, l_pos, l_nor, w, pdf):
    r_l, r_l_pos, r_l_nor, r_w_sum, r_M, r_pdf = r
    r_M += 1
    r_w_sum += w
    if random.random() * r_w_sum < w:
        r_l = l
        r_l_pos = l_pos
        r_l_nor = l_nor
        r_pdf = pdf
    return (r_l, r_l_pos, r_l_nor, r_w_sum, r_M, r_pdf)

class ReSTIR_DI_Engine(bpy.types.RenderEngine):
    # These three members are used by Blender to set up the
    # RenderEngine; define its internal name, visible name and capabilities.
    bl_idname = "ReSTIR_DI"
    bl_label = "ReSTIR_DI"
    bl_use_preview = False
    bl_use_shading_nodes = True
    bl_use_world_space_shading = True

    area_lights : List[AreaLight] = []
    ggx_mat_cache : GGXMaterialCache = None
    sampler : AreaLightsImportanceSampler = None
    reservoir = []
    ReSTIR_H = 10 # historical warm-up
    ReSTIR_M = 10 # new samples count
    ReSTIR_S_N = 10 # neightbors count
    ReSTIR_S_K = 5 # neightnors disk size
    ReSTIR_T_DECAY = 0.95
    ReSTIR_MAX_ENERGY = 50 # clamping weith to reduce fireflie

    # Init is called whenever a new render engine instance is created. Multiple
    # instances may exist at the same time, for example for a viewport and final
    # render.
    # Note the generic arguments signature, and the call to the parent class
    # `__init__` methods, which are required for Blender to create the underlying
    # `RenderEngine` data.
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scene_data = None
        self.draw_data = None

    # When the render engine instance is destroy, this is called. Clean up any
    # render engine data here, for example stopping running render threads.
    def __del__(self):
        # Own delete code...
        pass
    
    def update(self, data, depsgraph):
        # Minimal required for Blender to think the engine can evaluate materials
        pass
    
    def get_camera_details(self):
        cam = bpy.context.scene.camera
        cam_data = cam.data

        if cam_data.type == 'PERSP':
            if cam_data.sensor_fit == 'VERTICAL':
                sensor = cam.sensor_height
            else:  # HORIZONTAL or AUTO defaults to width
                sensor = cam_data.sensor_width
                
        cam_fov_rad = 2 * math.atan(sensor / (2 * cam_data.lens))
        cam_fov_deg = math.degrees(cam_fov_rad)
        cam_world_matrix = cam.matrix_world
        cam_pos = cam_world_matrix.translation
        cam_forward = cam_world_matrix.to_quaternion() @ Vector((0, 0, -1))
        cam_up = cam_world_matrix.to_quaternion() @ Vector((0, 1, 0))
        cam_right = cam_forward.cross(cam_up).normalized()
        
        ret = {
            "cam_fov_rad": cam_fov_rad,
            "cam_fov_deg": cam_fov_deg,
            "cam_world_matrix": cam_world_matrix,
            "cam_pos": cam_pos,
            "cam_forward": cam_forward,
            "cam_up": cam_up,
            "cam_right": cam_right
            }    
        
        return ret
        
    def generate_scene (self):
        def is_area_light(obj, ggx_mat_cache):
            for slot in obj.material_slots:
                k_d, k_s, k_r, k_e, k_es = ggx_mat_cache.get(slot.material)
                if k_es != 0 : return True
            return False
        self.area_lights = []
        self.ggx_mat_cache = GGXMaterialCache()
        depsgraph = bpy.context.evaluated_depsgraph_get()
        for instance in depsgraph.object_instances:
            obj = instance.object
            type = obj.type
            if type != 'MESH': continue
            if is_area_light(obj, self.ggx_mat_cache):
                self.area_lights.append(AreaLight(instance, self.ggx_mat_cache))

    def ray_trace (self, s, wi, d):
        scene = bpy.context.scene
        depsgraph = bpy.context.evaluated_depsgraph_get()
        result, pos, nor, index, object, matrix = scene.ray_cast(depsgraph, s, wi, distance = d)
        if result:
            mesh = object.data
            mat = self.ggx_mat_cache.get(object.material_slots[mesh.polygons[index].material_index])
            return result, pos, nor, mat
        return False, None, None, None
    
    def segment_visibility (self, s, e):
        wi = e-s
        d = wi.length
        wi /= d
        hit, p, n, m = self.ray_trace(s, wi, d)
        if not hit: return True
        return (p-e).length_squared < SAMPLER_BIAS * SAMPLER_BIAS * 9
    
    def restir_di (self, x, y, p, n, m, wo):
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0: return k_e * k_es

        def bdrf(wi):
            k_d, k_s, k_r, k_e, k_es = m
            return lambert_BDRF(k_d) + ggx_BDRF (wi, wo, n, Vector((0.04, 0.04, 0.04)), k_r) * k_s
        
        def target(sample):
            wi, pos, nor, pdf, cos_theta, l = sample
            return bdrf(wi) * l.li * cos_theta

        def target_pdf_l_p_n(sample):
            wi, pos, nor, pdf, cos_theta, l = sample
            return bdrf(wi) * l.li * cos_theta, pdf, l, pos, nor
        
        def visible(sample):
            wi, pos, nor, pdf, cos_theta, l = sample
            return self.segment_visibility(p + n * SAMPLER_BIAS, pos)

        def get_sample(reservoir):
            l, l_pos, l_nor, w_sum, M, pdf = reservoir
            if M == 0: return None
            pdf, wi = l.pdf(p, n, l_pos, l_nor)
            cos_theta = n.dot(wi)
            return (wi, l_pos, l_nor, pdf, cos_theta, l)
        
        def init_reservoir():
            r = reservoir_new()
            S = self.sampler.samples (p, n, self.ReSTIR_M)
            for s in S:
                if not visible(s): continue

                t, pdf, l, pos, nor = target_pdf_l_p_n(s)
                w = t.length / max(pdf, 1e-6)
                r = reservoir_update(r, l, pos, nor, w, pdf)
            return r
        
        def merge_reservoir(r, o_x, o_y, alpha):
            r_l, r_l_pos, r_l_nor, r_w_sum, r_M, r_pdf = r
            o_l, o_l_pos, o_l_nor, o_w_sum, o_M, o_pdf = self.reservoir[o_y*self.size_x+o_x]
            if o_M == 0 or o_w_sum == 0.0: return r

            # Compute target from other reservoir
            wi = (o_l_pos - p).normalized()
            cos_theta = n.dot(wi)
            if cos_theta <= 0: return r
            li = o_l.li
            t = li * bdrf(wi) * cos_theta

            w = (t.length / max(o_pdf, 1e-6)) * o_M
            w *= alpha

            r_M += 1
            r_w_sum += w

            if random.random() * r_w_sum < w:
                r_l = o_l
                r_l_pos = o_l_pos
                r_l_nor = o_l_nor
                r_pdf = t.length

            return (r_l, r_l_pos, r_l_nor, r_w_sum, r_M, r_pdf)
        
        
        def shade(reservoir):
            l, l_pos, l_nor, w_sum, M, pdf = reservoir
            if M == 0.0 or w_sum == 0.0: return Vector ((0,0, 0))
            s = get_sample(reservoir)
            t = target(s)
            return (t / max(pdf, 1e-6)) * (w_sum / M)
        
        r = init_reservoir()
        r = merge_reservoir(r, x, y, self.ReSTIR_T_DECAY)
        for _ in range(self.ReSTIR_S_N):
            radius = self.ReSTIR_S_K * math.sqrt(random.uniform(0, 1))
            angle = random.uniform(0, 2*math.pi)
            nx = int(x + math.cos(angle) * radius)
            ny = int(y + math.sin(angle) * radius)
            if nx < 0 or nx >= self.size_x or ny < 0 or ny >= self.size_y: continue
            r = merge_reservoir(r, nx, ny, 1.0)
        self.reservoir[y * self.size_x + x] = r
        return shade(r)

    # This is the method called by Blender for both final renders (F12) and
    # small preview for materials, world and lights.
    def render(self, depsgraph):
        scene = depsgraph.scene
        scale = scene.render.resolution_percentage / 100.0
        w = self.size_x = int(scene.render.resolution_x * scale)
        h = self.size_y = int(scene.render.resolution_y * scale)

        color = [0.0, 0.0, 0.0, 1.0]
        image = [color] * self.size_x * self.size_y
        
        # get the camera details
        cam = self.get_camera_details()
        cam_fov_rad = cam['cam_fov_rad']
        cam_pos = cam['cam_pos']
        cam_forward = cam['cam_forward']
        cam_up = cam['cam_up']
        cam_right = cam['cam_right']
        ar = self.size_x / self.size_y
        hty = math.tan(cam_fov_rad / 2)
        htx = hty * ar

        self.generate_scene()
        self.sampler = AreaLightsImportanceSampler(self.area_lights)
        r_empty = reservoir_new()
        self.reservoir = [r_empty for _ in range(w*h)]

        start = perf_counter()
        last_display = start
        self.ReSTIR_H = 2
        self.ReSTIR_M = 4
        self.ReSTIR_T_DECAY = 0.95
        self.ReSTIR_S_N = 4
        self.ReSTIR_S_K = w * 0.0625
        self.ReSTIR_MAX_ENERGY = 50 # clamping weith to reduce fireflie
        print (f"Rendering (  )...")
        print (f" - {self.ReSTIR_M} new samples from {self.sampler.get_name()}")
        print (f" - temporal warmp up {self.ReSTIR_H}")
        print (f" - temporal join decay {self.ReSTIR_T_DECAY}")
        print (f" - spatial radius {self.ReSTIR_S_K} pixels with {self.ReSTIR_S_N} samples")
        print (f" - max energy to reduce fireflies {self.ReSTIR_MAX_ENERGY}")
        for i in range(self.ReSTIR_H):

            for y in range(h):
                for x in range (w):
                    px = (2 * (x + 0.5) / w - 1) * htx
                    py = (1 - 2 * (y + 0.5) / h) * hty
                    wi = (cam_forward + cam_right * px + cam_up * py).normalized()
                    hit, p, n, m = self.ray_trace(cam_pos, wi, SAMPLER_DISTANCE)
                    if hit:
                        l = self.restir_di(x, y, p, n, m, wi)
                        image[(h-1-y)*w+x] = [l.x, l.y, l.z, 1]
                    now = perf_counter()
                    if now - last_display > 4.0:
                        last_display = now
                        ratio = (i * w * h + (y * w + x)) / (w*h*self.ReSTIR_H)
                        elapsed = now-start
                        total = elapsed / ratio
                        print (f"\tRender {int(100*ratio)} % Elapsed {get_cosmetic_duration(elapsed)} total {get_cosmetic_duration(total)} ETA {get_cosmetic_duration(total - elapsed)}")
        print (f"Rendering done in {get_cosmetic_duration(perf_counter() - start)}")
        
        self.area_lights = []
        self.ggx_mat_cache = None
        
        # Here we write the pixel values to the RenderResult
        result = self.begin_result(0, 0, self.size_x, self.size_y)
        layer = result.layers[0].passes["Combined"]
        layer.rect = image
        self.end_result(result)

    # For viewport renders, this method gets called once at the start and
    # whenever the scene or 3D viewport changes. This method is where data
    # should be read from Blender in the same thread. Typically a render
    # thread will be started to do the work while keeping Blender responsive.
    def view_update(self, context, depsgraph):
        region = context.region
        view3d = context.space_data
        scene = depsgraph.scene

        # Get viewport dimensions
        dimensions = region.width, region.height

        if not self.scene_data:
            # First time initialization
            self.scene_data = []
            first_time = True

            # Loop over all datablocks used in the scene.
            for datablock in depsgraph.ids:
                pass
        else:
            first_time = False

            # Test which datablocks changed
            for update in depsgraph.updates:
                print("Datablock updated: ", update.id.name)

            # Test if any material was added, removed or changed.
            if depsgraph.id_type_updated('MATERIAL'):
                print("Materials updated")

        # Loop over all object instances in the scene.
        if first_time or depsgraph.id_type_updated('OBJECT'):
            for instance in depsgraph.object_instances:
                pass

    # For viewport renders, this method is called whenever Blender redraws
    # the 3D viewport. The renderer is expected to quickly draw the render
    # with OpenGL, and not perform other expensive work.
    # Blender will draw overlays for selection and editing on top of the
    # rendered image automatically.
    def view_draw(self, context, depsgraph):
        # Lazily import GPU module, so that the render engine works in
        # background mode where the GPU module can't be imported by default.
        import gpu

        region = context.region
        scene = depsgraph.scene

        # Get viewport dimensions
        dimensions = region.width, region.height

        # Bind shader that converts from scene linear to display space,
        gpu.state.blend_set('ALPHA_PREMULT')
        self.bind_display_space_shader(scene)

        if not self.draw_data or self.draw_data.dimensions != dimensions:
            self.draw_data = ReSTIR_DI_DrawData(dimensions)

        self.draw_data.draw()

        self.unbind_display_space_shader()
        gpu.state.blend_set('NONE')

class ReSTIR_DI_DrawData:
    def __init__(self, dimensions):
        import gpu

        # Generate dummy float image buffer.
        self.dimensions = dimensions
        width, height = dimensions

        pixels = width * height * array.array('f', [0.1, 0.2, 0.1, 1.0])
        pixels = gpu.types.Buffer('FLOAT', width * height * 4, pixels)

        # Generate texture.
        self.texture = gpu.types.GPUTexture((width, height), format='RGBA16F', data=pixels)

        # Note: This is just a didactic example.
        # In this case it would be more convenient to fill the texture with:
        # self.texture.clear('FLOAT', value=[0.1, 0.2, 0.1, 1.0])

    def __del__(self):
        del self.texture

    def draw(self):
        from gpu_extras.presets import draw_texture_2d
        draw_texture_2d(self.texture, (0, 0), self.texture.width, self.texture.height)


# RenderEngines also need to tell UI Panels that they are compatible with.
# We recommend to enable all panels marked as BLENDER_RENDER, and then
# exclude any panels that are replaced by custom panels registered by the
# render engine, or that are not supported.
def get_panels():
    exclude_panels = {
        'VIEWLAYER_PT_filter',
        'VIEWLAYER_PT_layer_passes',
    }

    panels = []
    for panel in bpy.types.Panel.__subclasses__():
        if hasattr(panel, 'COMPAT_ENGINES') and 'BLENDER_RENDER' in panel.COMPAT_ENGINES:
            if panel.__name__ not in exclude_panels:
                panels.append(panel)

    return panels

def register():
    # Register the RenderEngine.
    bpy.utils.register_class(ReSTIR_DI_Engine)

    for panel in get_panels():
        panel.COMPAT_ENGINES.add('MonteCarlo')

def unregister():
    bpy.utils.unregister_class(ReSTIR_DI_Engine)

    for panel in get_panels():
        if 'MonteCarlo' in panel.COMPAT_ENGINES:
            panel.COMPAT_ENGINES.remove('MonteCarlo')

if __name__ == "__main__":
    register()