import math
import time
import random

# ============
# Math helpers
# ============

def normalize(v):
    l = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return (v[0]/l, v[1]/l, v[2]/l)

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
    
def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def mul_s(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def mul(a, b):
    return (a[0]*b[0], a[1]*b[1], a[2]*b[2])

def length(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def length_sq(v):
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]

def cross(a, b):
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )

def interp_barycentric(n0, n1, n2, u, v):
    w = 1.0 - u - v

    nx = w*n0[0] + u*n1[0] + v*n2[0]
    ny = w*n0[1] + u*n1[1] + v*n2[1]
    nz = w*n0[2] + u*n1[2] + v*n2[2]

    return (nx, ny, nz)


# ===========
# Ray tracing
# ===========

def ray_triangle_intersect(orig, dir, v0, v1, v2):
    EPSILON = 1e-6
    edge1 = sub(v1,v0)
    edge2 = sub(v2,v0)
    h = cross(dir,edge2)
    a = dot(edge1,h)
    if -EPSILON < a < EPSILON:
        return False, None, None, None
    f = 1.0 / a
    s = sub(orig, v0)
    u = f * dot(s,h)
    if u < 0.0 or u > 1.0:
        return False, None, None, None
    q = cross(s, edge1)
    v = f * dot(dir, q)
    if v < 0.0 or u + v > 1.0:
        return False, None, None, None
    t = f * dot(edge2,q)
    if t > EPSILON:
        return True, t, u, v
    return False, None, None, None


# Visibility check

def ray_visibility_face(p, d, dist, verts, inds):
    num = len(inds)
    for i in range(num - 2):
        a = verts[inds[0]]
        b = verts[inds[i + 1]]
        c = verts[inds[i + 2]]
        hit, t, u, v = ray_triangle_intersect(p, d, a, b, c)
        if hit and t <= dist: return True
    return False

def ray_visibility_mesh(p, d, dist, mesh, ignore_emissive):
    (verts, mats, face_inds, face_nors, face_mats) = mesh

    for i, inds in enumerate(face_inds):
        nor = face_nors[i]
        if dot(nor, d) > 0.0: continue
        if ignore_emissive:
            mat = mats[face_mats[i]]
            c, ems = mat
            if length(ems) > 0.0: continue
        if ray_visibility_face(p, d, dist, verts, inds): return True

    return False

def ray_visibility_world(p, d, dist, meshes, ignore_emissive):
    for mesh in meshes:
        if ray_visibility_mesh(p, d, dist, mesh, ignore_emissive): return True
    return False

def segment_visibility_world(s, e, meshes, ignore_emissive):
    d = sub (e, s)
    dist = length(d)
    d = mul_s(d, 1.0 / dist)
    return ray_visibility_world(s, d, dist, meshes, ignore_emissive)

# trace

def ray_trace_face(p, d, verts, inds):
    num = len(inds)
    for i in range(num - 2):
        a = verts[inds[0]]
        b = verts[inds[i + 1]]
        c = verts[inds[i + 2]]
        hit, dist, u, v = ray_triangle_intersect(p, d, a, b, c)
        if hit: 
            return True, dist
    return False, None

def ray_trace_mesh(p, d, dist, mesh):
    (verts, mats, face_inds, face_nors, face_mats) = mesh

    hit_dist = dist
    hit_face_index = None
    hit_mat = None

    for i, inds in enumerate(face_inds):
        nor = face_nors[i]
        if dot(nor, d) > 0.0: continue
        hit, t = ray_trace_face(p, d, verts, inds)
        if hit and t < hit_dist:
            hit_dist = t
            hit_face_index = i
            hit_mat = mats[face_mats[i]]

    if hit_dist != dist:
        o = mul_s(d, hit_dist)
        pos = add(p, o)
        nor = face_nors[hit_face_index]
        return True, pos, nor, hit_dist, hit_face_index, hit_mat

    return False, None, None, None, None, None

def ray_trace_world(p, d, meshes):
    hit_dist = math.inf
    hit_pos = None
    hit_nor = None
    hit_mat = None
    for mesh in meshes:
        hit, pos, nor, t, face, mat = ray_trace_mesh(p, d, hit_dist, mesh)
        if hit and t < hit_dist:
            hit_dist = t
            hit_pos = pos
            hit_nor = nor
            hit_mat = mat
    return hit_dist != math.inf, hit_pos, hit_nor, hit_mat

def ray_trace_world_ex(p, d, meshes):
    hit_dist = math.inf
    hit_pos = None
    hit_nor = None
    hit_mat = None
    hit_mesh_index = None
    hit_face_index = None
    for mesh_index, mesh in enumerate(meshes):
        hit, pos, nor, t, face, mat = ray_trace_mesh(p, d, hit_dist, mesh)
        if hit and t < hit_dist:
            hit_dist = t
            hit_pos = pos
            hit_nor = nor
            hit_mat = mat
            hit_mesh_index = mesh_index
            hit_face_index = face
    return hit_dist != math.inf, hit_pos, hit_nor, hit_mat, hit_mesh_index, hit_face_index




# ========
# samplers
# ========

def get_orthonormal_axis(z_axis):
    # generate perp vectors
    x_axis = (1.0, 0.0, 0.0)
    y_axis = cross(z_axis, x_axis)
    if length(y_axis) < 0.01:
        x_axis = (0.0, 1.0, 0.0)
        y_axis = cross(z_axis, x_axis)
        
    x_axis = cross(y_axis, z_axis)

    return x_axis, y_axis

def get_area_light(area_lights, mesh_index, face_index):
    for light in area_lights:
        light_mesh_index, area, ems_weight, faces = light
        if mesh_index != light_mesh_index: continue
        for face in faces:
            area, ems_weight, light_face_index, areas = face
            if face_index != light_face_index: continue
            return light
    return None

class UniformSampler:
    meshes = None
    area_lights = None

    def __init__(self, meshes, area_lights):
        self.meshes = meshes
        self.area_lights = area_lights

    def sample (self, p, n):
        x_axis, y_axis = get_orthonormal_axis(n)
        u = random.uniform(0.0, 2.0 * math.pi)
        cos_v = random.uniform(0.0, 1.0)
        v = math.acos(cos_v)

        x = math.sin(v) * math.cos(u)
        y = math.sin(v) * math.sin(u) 
        z = math.cos(v)

        wi = mul_s(x_axis, x)
        wi = add(wi, mul_s(y_axis, y))
        wi = add(wi, mul_s(n, z))

        hit, pos, nor, mat, mesh_index, face_index = ray_trace_world_ex (p, wi, self.meshes)
        area = None
        weight = None
        li = (0, 0, 0)
        if hit: 
            _, li = mat
            if length(li) > 0.0:
                _, area, weight, _ = get_area_light(self.area_lights, mesh_index, face_index)
            
        
        pdf = 1.0 / (2.0 * math.pi)

        cos_theta = z

        return li, pdf, cos_theta, pos, nor, area, weight

class CosineSampler:
    meshes = None
    area_lights = None

    def __init__(self, meshes, area_lights):
        self.meshes = meshes
        self.area_lights = area_lights

    def sample (self, p, n):
        x_axis, y_axis = get_orthonormal_axis(n)
        u = random.uniform(0.0, 1.0)
        v = random.uniform(0.0, 2.0 * math.pi)

        r = math.sqrt(u)
        x = r*math.cos(v)
        y = r*math.sin(v)
        z = math.sqrt(1.0-u)

        wi = mul_s(x_axis, x)
        wi = add(wi, mul_s(y_axis, y))
        wi = add(wi, mul_s(n, z))

        hit, pos, nor, mat, mesh_index, face_index = ray_trace_world_ex (p, wi, self.meshes)
        area = None
        weight = None
        li = (0, 0, 0)
        if hit: 
            _, li = mat
            if length(li) > 0.0:
                _, area, weight, _ = get_area_light(self.area_lights, mesh_index, face_index)

        pdf = z / math.pi

        cos_theta = z

        return li, pdf, cos_theta, pos, nor, area, weight
    
class AreaLightsSampler:
    meshes = None
    area_lights = None
    def __init__(self, meshes, area_lights):
        self.meshes = meshes
        self.area_lights = area_lights

    def __get_weighted_random_item (items, total_weight, get_item_weight, use_index = False):
        weight = random.random() * total_weight
        accumulator = 0.0
        for index, i in enumerate(items):
            accumulator = accumulator + get_item_weight(i)
            if accumulator >= weight: 
                if use_index: return index
                return i
        if use_index: return len(items)-1
        return items[-1]
    
    def compute_pdf(p, pos, nor, area, weight, total_weight):
        pdf_area_uniform = 1.0 / area
        wi = sub(pos,p)
        r = length(wi)
        wi = mul_s(wi, 1.0 / r)
        cos_theta_light =  -dot(nor, wi)
        if cos_theta_light <= 0.0: return 0.0
        pdf = pdf_area_uniform * r * r / cos_theta_light
        pdf *= weight / total_weight
        return pdf

    def sample(self,p,n):
        # get random sample on weight random light
        total_ems_weight = self.area_lights[-1]
        mesh_index, area, ems_weight, faces = AreaLightsSampler.__get_weighted_random_item(self.area_lights, total_ems_weight, lambda light: light[2])
        face_area, face_ems_weight, face_index, face_areas = AreaLightsSampler.__get_weighted_random_item(faces, ems_weight, lambda face: face[1])
        triangle = AreaLightsSampler.__get_weighted_random_item(face_areas, face_area, lambda t:t, use_index=True)

        verts, mats, face_inds, face_nors, face_mats = self.meshes[mesh_index]
        inds = face_inds[face_index]

        a = verts[inds[0]]
        b = verts[inds[triangle + 1]]
        c = verts[inds[triangle + 2]]
        u = random.random()
        v = random.random()
        while u + v > 1.0: v = random.random()

        sample = interp_barycentric(a, b, c, u, v)
        nor = face_nors[face_index]
        _, li = mats[face_mats[face_index]]
        wi = sub(sample,p)
        r = length(wi)
        wi = mul_s(wi, 1.0 / r)

        cos_theta = max(0, dot(n, wi))

        # check visibility
        pos = sample
        if segment_visibility_world(p, sample, self.meshes, True):  
            li = (0, 0, 0)

        # compute the pdf light sample
        pdf_area_uniform = 1.0 / area # pdf for sample in respect to area, we need solid angle pdf so
        cos_theta_light = -dot(nor, wi)
        pdf = pdf_area_uniform * r * r / cos_theta_light

        # adjust the light sample pdf with the pdf of the light
        pdf *= ems_weight / total_ems_weight

        return li, pdf, cos_theta, pos, nor, area, ems_weight
    
g_uniform_sampler = None
g_cosine_sampler = None
g_area_light_sampler = None

def init_samplers(meshes, area_lights):
    global g_uniform_sampler, g_cosine_sampler, g_area_light_sampler
    g_uniform_sampler = UniformSampler(meshes, area_lights)
    g_cosine_sampler = CosineSampler(meshes, area_lights)
    g_area_light_sampler = AreaLightsSampler(meshes, area_lights)

def destroy_samplers():
    global g_uniform_sampler, g_cosine_sampler, g_area_light_sampler
    del g_uniform_sampler
    del g_cosine_sampler
    del g_area_light_sampler



# ===========================================================
# Radiance via Monte Carlor solving of the rendering equation
# ===========================================================

# radiance output from solving the rendering equation with 1 Monte Carlo Sampler and Lamber BDRF
def radiance_from_sampler(num_samples, p, n, sampler):
    accumulator = (0, 0, 0)
    for i in range(num_samples):
        li, pdf, cos_theta, pos, nor, area, weight = sampler.sample(p, n)
        
        brdf = 1.0 / math.pi # lambert

        accumulator = add(accumulator, mul_s(li, brdf * cos_theta / pdf))

    return mul_s(accumulator, 1.0 / num_samples)

# MIS
def mis_weight(pdf_i, pdf_1, pdf_2, n=2):
    return (pdf_i ** n) / (pdf_1 ** n + pdf_2 ** n)

def radiance_from_mis(num_samples, p, n, sampler1, sampler2):
    accumulator = (0, 0, 0)
    brdf = 1.0 / math.pi # lambert
    for i in range(num_samples):

        li1, pdf1, cos_theta1, p1, n1, area1, weight1 = sampler1.sample(p, n)
        li2, pdf2, cos_theta2, p2, n1, area2, weight2 = sampler2.sample(p, n)
        
        if cos_theta1 > 0:
            w = mis_weight(pdf1, pdf1, pdf2)
            accumulator = add(accumulator, mul_s(li1, brdf * cos_theta1 * w  / pdf1))

        if cos_theta2 > 0:
            w = mis_weight(pdf2, pdf2, pdf1)
            accumulator = add(accumulator, mul_s(li2, brdf * cos_theta2 * w  / pdf2))

    return mul_s(accumulator, 1.0 / num_samples)

# =======================
# Multi process rendering
# =======================
    
# Generate the image pixel by pixel
def compute_pixel_color(x, y, tracing_context):
    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context
    
    # generate the ray dir
    px = (2 * (x + 0.5) / w - 1) * htx
    py = (1 - 2 * (y + 0.5) / h) * hty
    h = mul_s(cam_right, px)
    v = mul_s(cam_up, py)
    d = add(cam_forward, h)
    d = add(d, v)
    d = normalize(d)

    color = (0.0, 0.0, 0.0)

    # search for the point
    hit, pos, nor, mat = ray_trace_world (cam_pos, d, meshes)
    if not hit: return [color[0], color[1], color[2], 1.0]
    (c, li) = mat
    if length(li) > 0.0: return [li[0], li[1], li[2], 1.0]

    # get the incoming light
    l = radiance_from_mis(1, pos, nor, g_cosine_sampler, g_area_light_sampler)
#    l = radiance_from_sampler(20, pos, nor, g_cosine_sampler)
    color = mul(l, c)

    return [color[0], color[1], color[2], 1.0]


g_progress_lock = None
g_progress = None
g_progress_last_display = None

def init(progress_lock, progress, progress_last_display):
    global g_progress_lock, g_progress, g_progress_last_display
    g_progress_lock = progress_lock
    g_progress = progress
    g_progress_last_display = progress_last_display
    
def compute_process(args):
    x1, y1, x2, y2, tracing_context = args

    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context

    init_samplers(meshes, area_lights)

    local_buffer = []
    start = time.time()
    processed = 0
    for y in range(y1, y2):
        for x in range(x1, x2):
            local_buffer.append(compute_pixel_color(x, y, tracing_context))

            processed += 1
            now = time.time()
            if now - start > 2.0:
                start = now
                with g_progress_lock:
                    g_progress.value += processed
                    processed = 0
                    if now - g_progress_last_display.value > 2.0:
                        p = g_progress.value / (w * h)
                        print (f"\trendering {int(p * 100)}%")
                        g_progress_last_display.value = now

    destroy_samplers()
            
    return local_buffer

# ======
# ReSTIR
# ======

# a simplified restir lite method.
ReSTIR_TEMPORAL_NUM_SAMPLES = 1
ReSTIR_TEMPORAL_DROP_RATIO = 0.0
ReSTIR_SPATIAL_KERNEL = 7
ReSTIR_SPATIAL_DROP_RATIO = 0.75
ReSTIR_CHECK_VISBILITY_RATIO = 0.0


def get_empty_reservoir():
    p = None
    n = None
    mat = None
    samples = None
    return (p, n, mat, samples)

def generate_restir_reservoir (w, h):
    empty = get_empty_reservoir()
    return [empty for _ in range(h * w)]

def reservoir_add(reservoir, p, n, mat, samples): reservoir.append((p, n, mat, samples))

def init_pixel_reservoir(x, y, tracing_context, reservoir):
    global ReSTIR_TEMPORAL_NUM_SAMPLES
    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context
    
    # generate the ray dir
    px = (2 * (x + 0.5) / w - 1) * htx
    py = (1 - 2 * (y + 0.5) / h) * hty
    h = mul_s(cam_right, px)
    v = mul_s(cam_up, py)
    d = add(cam_forward, h)
    d = add(d, v)
    d = normalize(d)

    empty = get_empty_reservoir()

    hit, pos, nor, mat = ray_trace_world (cam_pos, d, meshes)
    if not hit > 0.0: 
        reservoir.append(empty)
        return
    
    c, li = mat
    if length(li) > 0.0:
        reservoir_add (reservoir, pos, nor, mat, []) 
    else:
        # get the incoming light
        reservoir_add (reservoir, pos, nor, mat, [g_area_light_sampler.sample(pos, nor) for _ in range(ReSTIR_TEMPORAL_NUM_SAMPLES)])

def compute_process_restir_init(args):
    x1, y1, x2, y2, tracing_context = args

    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context

    init_samplers(meshes, area_lights)

    local_reservoir = []
    start = time.time()
    processed = 0
    for y in range(y1, y2):
        for x in range(x1, x2):
            init_pixel_reservoir(x, y, tracing_context, local_reservoir)

            processed += 1
            now = time.time()
            if now - start > 2.0:
                start = now
                with g_progress_lock:
                    g_progress.value += processed
                    processed = 0
                    if now - g_progress_last_display.value > 2.0:
                        p = g_progress.value / (w * h)
                        print (f"\tInit ReSTIR {int(p * 100)}%")
                        g_progress_last_display.value = now

    destroy_samplers()
            
    return local_reservoir


def neighbor_match(p, n, np, nn):
    if not p or not np: return False
    if dot(n, nn) < 0.6: return False
    if length_sq(sub(p, np)) > 0.5: return False
    return True


def get_pixel_variance(x,y, w, h, kernel, reservoir):
    values = []
    p, n, mat, samples = reservoir[y * w + x]
    for i in range(y-kernel, y+kernel+1):
        for j in range(x-kernel, x+kernel+1):
            if i < 0 or i >= h or j < 0 or j >= w: continue
            np, nn, mat, samples = reservoir[i * w + j]
            if not neighbor_match(p, n, np, nn) or length (mat[1]) > 0.0: continue

            for li, pdf, cos_theta, pos, nor, area, weight in samples:
                values.append (length(li) * cos_theta / pdf)
    
    if len(values) < 2: return 0.0

    mean = sum(values) / len(values)
    return sum((v - mean) ** 2 for v in values) / (len(values) - 1)

def compute_pixel_color_restir(x, y, tracing_context, reservoir):
    global ReSTIR_SPATIAL_KERNEL, ReSTIR_CHECK_VISBILITY_RATIO, ReSTIR_SPATIAL_DROP_RATIO, ReSTIR_TEMPORAL_DROP_RATIO
    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context
    total_weight = area_lights[-1]

    color = (0, 0, 0)
    p, n, mat, samples = reservoir[y * w + x]
    if not p: return [color[0], color[1], color[2], 1.0]

    c, li = mat
    if length(li) > 0:  return [li[0], li[1], li[2], 1.0]

    brdf = 1.0 / math.pi # lambert
    accumulator = (0, 0, 0)
    num_samples = len(samples)

    for li, pdf, cos_theta, pos, nor, area, weight in samples:
        accumulator = add(accumulator, mul_s(li, brdf * cos_theta / pdf))

    def get_kernel(x, y, w, h):
        variance = get_pixel_variance(x, y, w, h, 3, reservoir)
        if variance < 1.0: return ReSTIR_SPATIAL_KERNEL
        if variance < 10.0: return ReSTIR_SPATIAL_KERNEL // 2
        if variance < 500.0: return ReSTIR_SPATIAL_KERNEL // 4
        return ReSTIR_SPATIAL_KERNEL // 8

    for i in range(y-ReSTIR_SPATIAL_KERNEL, y+ReSTIR_SPATIAL_KERNEL+1):
        for j in range(x-ReSTIR_SPATIAL_KERNEL, x+ReSTIR_SPATIAL_KERNEL+1):
            if i < 0 or i >= h or j < 0 or j >= w: continue
            if i == y and j == x: continue
            if random.random() < ReSTIR_SPATIAL_DROP_RATIO: continue

            p2, n2, mat2, samples2 = reservoir[i * w + j]
            if not p2: continue

            if dot(n, n2) < 0.6: continue

            if length_sq(sub(p, p2)) > 0.5: continue

            c2, li2 = mat
            if length(li2) > 0:continue

            num_samples2 = len(samples2)
            for li, pdf, cos_theta, pos, nor, area, weight in samples2:
                if num_samples2 > 1 and random.random() < ReSTIR_TEMPORAL_DROP_RATIO: continue
                if not pos: continue

                num_samples += 1

                # recompute cos_theta
                wi = normalize(sub(pos, p))
                cos_theta = dot(wi, n)
                if cos_theta <= 0: continue

                # recompute pdf
                pdf = AreaLightsSampler.compute_pdf(p, pos, nor, area, weight, total_weight)
                if pdf <= 0.0: continue

                check_visibility = random.random() > 1.0 - ReSTIR_CHECK_VISBILITY_RATIO
                if check_visibility and segment_visibility_world(p, pos, meshes, True): continue

                accumulator = add(accumulator, mul_s(li, brdf * cos_theta / pdf))

    accumulator = mul_s(accumulator, 1.0 / num_samples)

    color = mul(c, accumulator)

    return [color[0], color[1], color[2], 1.0]

def compute_process_restir_render(args):
    x1, y1, x2, y2, tracing_context, reservoir = args

    (w, h, cam_pos, cam_forward, cam_right, cam_up, htx, hty, meshes, area_lights) = tracing_context

    init_samplers(meshes, area_lights)

    local_buffer = []
    start = time.time()
    processed = 0
    for y in range(y1, y2):
        for x in range(x1, x2):
            local_buffer.append(compute_pixel_color_restir(x, y, tracing_context, reservoir))

            processed += 1
            now = time.time()
            if now - start > 2.0:
                start = now
                with g_progress_lock:
                    g_progress.value += processed
                    processed = 0
                    if now - g_progress_last_display.value > 2.0:
                        p = g_progress.value / (w * h)
                        print (f"\tReSTIR rendering {int(p * 100)}%")
                        g_progress_last_display.value = now

    destroy_samplers()
            
    return local_buffer
