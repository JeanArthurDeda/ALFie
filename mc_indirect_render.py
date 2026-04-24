import sys
import os
import bpy
import array
import math
import random
from abc import ABC, abstractmethod
from mathutils import Vector, Matrix
from time import perf_counter
import ctypes
from typing import List

# Add current .blend directory to Python path
blend_dir = os.path.dirname(bpy.data.filepath)
if blend_dir not in sys.path:
    sys.path.append(blend_dir)
#from restir_helpers import generate_scene


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
    mat = None
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
        self.mat = self.ggx_mat_cache.get(obj.material_slots[0].material)
        k_d, k_s, k_r, k_e, k_es = self.mat

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

# GGX
def build_tangent_frame(n: Vector):
    if abs(n.z) < 0.999:
        t = Vector((0.0, 0.0, 1.0)).cross(n).normalized()
    else:
        t = Vector((1.0, 0.0, 0.0)).cross(n).normalized()
    b = n.cross(t)
    return t, b
def ggx_sample_vndf(wo: Vector, n: Vector, roughness: float):
    """
    GGX VNDF sampling (Heitz 2018).
    Returns: wi, pdf
    """

    # Transform wo to local space
    t, b = build_tangent_frame(n)
    wo_local = Vector((
        wo.dot(t),
        wo.dot(b),
        wo.dot(n)
    )).normalized()

    # Handle degenerate case
    if wo_local.z <= 0.0:
        return Vector((0, 0, 0)), 0.0, 0.0

    # GGX alpha
    alpha = roughness * roughness

    # Stretch view direction
    Vh = Vector((
        alpha * wo_local.x,
        alpha * wo_local.y,
        wo_local.z
    )).normalized()

    # Orthonormal basis
    lensq = Vh.x * Vh.x + Vh.y * Vh.y
    if lensq > 0.0:
        T1 = Vector((-Vh.y, Vh.x, 0.0)) / math.sqrt(lensq)
        T2 = Vh.cross(T1)
    else:
        T1 = Vector((1.0, 0.0, 0.0))
        T2 = Vector((0.0, 1.0, 0.0))

    # Sample point on disk
    u1 = random.random()
    u2 = random.random()

    r = math.sqrt(u1)
    phi = 2.0 * math.pi * u2

    t1 = r * math.cos(phi)
    t2 = r * math.sin(phi)

    # Warp t2
    s = 0.5 * (1.0 + Vh.z)
    t2 = (1.0 - s) * math.sqrt(max(0.0, 1.0 - t1 * t1)) + s * t2

    # Reproject onto hemisphere
    Nh = (t1 * T1 + t2 * T2 + math.sqrt(max(0.0, 1.0 - t1*t1 - t2*t2)) * Vh)

    # Unstretch
    h_local = Vector((
        alpha * Nh.x,
        alpha * Nh.y,
        max(0.0, Nh.z)
    )).normalized()

    # Transform back to world
    h = (t * h_local.x + b * h_local.y + n * h_local.z).normalized()

    # Reflect wo around h
    wi = 2.0 * wo.dot(h) * h - wo
    wi.normalize()

    # Reject below surface
    n_dot_wi = n.dot(wi)
    if n_dot_wi <= 0.0:
        return wi, 0.0, 0.0

    # ---- PDF (VNDF-consistent) ----
    n_dot_h = max(0.0, n.dot(h))
    wo_dot_h = max(0.0, wo.dot(h))
    n_dot_wo = max(0.0, n.dot(wo))

    if wo_dot_h <= 0.0 or n_dot_wo <= 0.0:
        return wi, 0.0, 0.0

    # GGX D (consistent alpha usage)
    alpha2 = alpha * alpha
    denom = (n_dot_h * n_dot_h) * (alpha2 - 1.0) + 1.0
    D = alpha2 / (math.pi * denom * denom)

    # Smith G1 (only for wo)
    def G1(n_dot_v):
        if n_dot_v <= 0.0:
            return 0.0
        tan2 = (1.0 - n_dot_v * n_dot_v) / (n_dot_v * n_dot_v)
        return 2.0 / (1.0 + math.sqrt(1.0 + alpha2 * tan2))

    G1_wo = G1(n_dot_wo)

    pdf = (D * G1_wo * wo_dot_h) / n_dot_wo

    return wi, pdf, n_dot_wi


def ggx_pdf_vndf(wi: Vector, wo: Vector, n: Vector, roughness: float) -> float:
    """
    PDF for GGX VNDF sampling (Heitz).

    Args:
        wi: incoming light direction (toward light)
        wo: outgoing/view direction (toward camera)
        n: surface normal
        roughness: perceptual roughness

    Returns:
        pdf value (solid angle measure)
    """

    wi = wi.normalized()
    wo = wo.normalized()
    n = n.normalized()

    n_dot_wi = n.dot(wi)
    n_dot_wo = n.dot(wo)

    # Must be above surface
    if n_dot_wi <= 0.0 or n_dot_wo <= 0.0:
        return 0.0

    # Half-vector
    h = (wi + wo)
    if h.length_squared == 0.0:
        return 0.0
    h.normalize()

    n_dot_h = max(0.0, n.dot(h))
    wo_dot_h = max(0.0, wo.dot(h))

    if n_dot_h <= 0.0 or wo_dot_h <= 0.0:
        return 0.0

    # --- GGX parameters ---
    alpha = roughness * roughness
    alpha2 = alpha * alpha

    # --- GGX NDF (D) ---
    denom = (n_dot_h * n_dot_h) * (alpha2 - 1.0) + 1.0
    D = alpha2 / (math.pi * denom * denom)

    # --- Smith G1 for wo only ---
    def G1(n_dot_v):
        if n_dot_v <= 0.0:
            return 0.0
        tan2 = (1.0 - n_dot_v * n_dot_v) / (n_dot_v * n_dot_v)
        return 2.0 / (1.0 + math.sqrt(1.0 + alpha2 * tan2))

    G1_wo = G1(n_dot_wo)

    # --- VNDF PDF ---
    pdf = (D * G1_wo * wo_dot_h) / n_dot_wo

    return max(0.0, pdf)

SAMPLER_DISTANCE = 20.0
SAMPLER_BIAS = 0.001
class Sampler(ABC):
    # in case the sampler supports sampling from a unnormalized uniform pdf or target function. RIS supports sampling from target pdf
    # target_pdf (wi, cos_theta) -> float
    param_target_pdf = None
    # some samplers (GGX) needs the current material
    param_mat = None

    samplers = []

    def set_params(self, param_target_pdf, param_mat):
        self.param_target_pdf = param_target_pdf
        self.param_mat = param_mat
        for s in self.samplers:
            s.set_params(param_target_pdf, param_mat)

    # return [(wi, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    @abstractmethod
    def samples(self, p, n, wo, num):
        pass
    # return pdf
    @abstractmethod
    def pdf(self, p, n, wo, wi):
        pass
    @abstractmethod
    def get_name(self):
        pass

class SimpleSampler(Sampler):
    sample_l = None
    pdf_l = None
    name = None

    def __init__(self, sample_l, pdf_l, name):
        self.sample_l = sample_l
        self.pdf_l = pdf_l
        self.name = name

    def samples (self, p, n, wo, num):
        S = []
        for _ in range(num):
            wi, pdf, cos_theta = self.sample_l(n)
            S.append((wi, pdf, 1.0, cos_theta))
        return S
    
    def pdf(self, p, n, wo, wi):
        return self.pdf_l(p, n, wi)

    def get_name(self):
        return self.name

class AreaLightsImportanceSampler(Sampler):
    area_lights : List[AreaLight] = []
    weight = 0

    def __init__(self, area_lights : List[AreaLight]):

        self.area_lights = area_lights
        for l in self.area_lights:
            self.weight += l.weight

    def samples(self, p, n, wo, num):
        S = []
        for _ in range(num):
            l = self.area_lights[get_weighted_random_index(self.area_lights, self.weight, lambda l: l.weight)]
            pos, nor = l.get_weight_random_point()
            pdf, wi = l.pdf(p, n, pos, nor)
            pdf *= l.weight / self.weight # adjust pdf based on light selection pdf
            cos_theta = max(0, wi.dot(n))
            S.append((wi, pdf, 1.0, cos_theta))
        return S
    
    def pdf(self, p, n, wo, wi):
        pdf = 0.0
        for l in self.area_lights:
            hit, pos, nor = l.ray_trace(p, wi, SAMPLER_DISTANCE)
            if not hit: continue
            l_pdf, _ = l.pdf(p, n, pos, nor)
            pdf += l_pdf * l.weight / self.weight
        return pdf
    
    def get_name(self):
        return f"LightsImportance({len(self.area_lights)} lights)"
    
class GGXSampler(Sampler):
    def __init__(self):
        super().__init__()

    def samples(self, p, n, wo, num):
        S = []
        k_d, k_s, k_r, k_e, k_es = self.param_mat
        for _ in range(num):
            wi, pdf, cos_theta = ggx_sample_vndf(wo, n, k_r)
            S.append ((wi, pdf, 1.0, cos_theta))
        return S
    
    def pdf(self, p, n, wo, wi):
        if not self.param_mat: return 0.0
        k_d, k_s, k_r, k_e, k_es = self.param_mat
        return ggx_pdf_vndf(wi, wo, n, k_r)
    
    def get_name(self):
        return "GGX"
    
class MISSampler(Sampler):
    s1 : Sampler = None
    s2 : Sampler = None
    ratio : float = None

    def __init__(self, s1, s2, ratio):
        super().__init__()
        self.samplers = [s1, s2]
        self.ratio = ratio

    def samples(self, p, n, wo, num):
        r = num * self.ratio
        num1 = int(r)
        if random.random() < (r - num1): num1 += 1
        num2 = num - num1
        s1 = self.samplers[0]
        s2 = self.samplers[1]
        S1 = s1.samples(p,n, wo, num1)
        S2 = s2.samples(p,n, wo, num2)
        S = []
        for s in S1:
            wi1, pdf1, mis_w1, cos_theta1 = s
            if pdf1 == 0: 
                S.append(s)
                continue
            pdf2 = s2.pdf(p, n, wo, wi1)
            w = (num1 * pdf1) / (num1 * pdf1 + num2 * pdf2)
            S.append((wi1, pdf1, w, cos_theta1))
        for s in S2:
            wi2, pdf2, mis_w2, cos_theta2 = s
            if pdf2 == 0: 
                S.append(s)
                continue
            pdf1 = s1.pdf(p, n, wo, wi2)
            w = (num2 * pdf2) / (num1 * pdf1 + num2 * pdf2)
            S.append((wi2, pdf2, w, cos_theta2))
        return S
    
    def pdf(self, p, n, wo, wi): # MIS Samples cannot be used in any strategy that requires computing the PDF for a generic wi
        return None
    
    def get_name(self):
        s1 = self.samplers[0]
        s2 = self.samplers[1]
        return f"MIS {int(self.ratio * 100)}% {s1.get_name()} {100 - int(self.ratio * 100)}% {s2.get_name()}"
    
class RISSampler(Sampler):
    M : int = None

    def __init__(self, M : int, s : Sampler):
        super().__init__()
        self.M = M
        self.samplers = [s]

    def samples(self, p, n, wo, num):
        F = []

        s = self.samplers[0]
        for _ in range(num):
            S = s.samples(p, n, wo, self.M)
            W = [0.0 if pdf == 0.0 or mis_w == 0.0 else self.param_target_pdf(wi, cos_theta) * mis_w / pdf for wi, pdf, mis_w, cos_theta in S]
            total_weight = sum(W)

            champion = get_weighted_random_index(W, total_weight, lambda w: w)
            wi, pdf, mis_w, cos_theta = S[champion]
            if W[champion] == 0.0:
                F.append((wi, 0.0, 0.0, cos_theta))
            else:
                f = W[champion] * pdf / mis_w # compute target_pdf from weight instead of calculating it
                w = (total_weight / self.M) / f
                # RIS replaces the mis_w as 1 it's already incorporated into w. We can have RIS from MIS
                # and 1/w can be used as a PDF in the sense that f(x) / pdf = f(x) * w
                F.append((wi, 1/w, 1.0, cos_theta)) 

        return F
    
    def pdf(self, p, n, wo, wi):  # MIS Samples cannot be used in any strategy that requires computing the PDF for a generic wi
        return None
    
    def get_name(self):
        s = self.samplers[0]
        return f"RIS({self.M} from {s.get_name()})"

        

# BDRF =========================

def lambert_BDRF(k_d):
    return k_d / math.pi

def phong_BRDF(wi, wo, N):
    k_s = 0.75
    shininess = 32.0
    
    R = 2.0 * N.dot(wi) * N - wi
    
    specular = k_s * math.pow(max(0, R.dot(wo)), shininess) #(shininess + 2)/(2*math.pi)
    
    return specular

def ggx_BDRF(wi: Vector, wo: Vector, N: Vector, F0: Vector, roughness: float) -> float:
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


class MonterCarloIndirectEngine(bpy.types.RenderEngine):
    # These three members are used by Blender to set up the
    # RenderEngine; define its internal name, visible name and capabilities.
    bl_idname = "MonteCarloGI"
    bl_label = "MonteCarloGI"
    bl_use_preview = False
    bl_use_shading_nodes = True
    bl_use_world_space_shading = True

    area_lights : List[AreaLight] = []
    ggx_mat_cache : GGXMaterialCache = None
    direct_sampler : Sampler = None
    indirect_sampler : Sampler = None
    max_bounce = 1
    russia_roulette_level = 1
    russian_roulette_prob = 0.5
    debug_points = []

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
    
    def compute_randiance (self, p, n, m, wo, num, level):
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0 or level > self.max_bounce:
            return k_e * k_es
        
        # russian roulette
        russian_roulette_scale = 1.0
        if level >= self.russia_roulette_level:
            if random.random() > self.russian_roulette_prob: return k_e * k_es
            russian_roulette_scale = 1.0 / self.russian_roulette_prob

        def bdrf(wi):
            k_d, k_s, k_r, k_e, k_es = m
            return lambert_BDRF(k_d) + ggx_BDRF (wi, wo, n, Vector((0.04, 0.04, 0.04)), k_r) * k_s

        accumulator = Vector((0, 0, 0))
        sampler = self.direct_sampler if level == 0 else self.indirect_sampler
        sampler.set_params( lambda wi, cos_theta: bdrf (wi).length * cos_theta, m)
        S = sampler.samples (p, n, wo, num)
        for s in S:
            wi, pdf, mis_w, cos_theta = s
            if pdf == 0.0: continue
            hit, pos, nor, mat = self.ray_trace(p + n * SAMPLER_BIAS, wi, SAMPLER_DISTANCE)
            if not hit: continue
            li = self.compute_randiance(pos, nor, mat, -wi, num, level + 1)
            accumulator += bdrf(wi) * li * mis_w * cos_theta / pdf
        accumulator = russian_roulette_scale * accumulator / num
        return accumulator
    
    def primary (self, p, wi, num):
        hit, pos, nor, mat = self.ray_trace(p, wi, SAMPLER_DISTANCE)
        if hit:
            l = self.compute_randiance(pos, nor, mat, -wi, num, 0)
            return [l.x, l.y, l.z, 1]
        return [0, 0, 0, 0]

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
        cosine = SimpleSampler(cosine_sample, cosine_pdf, "Cosine")
        area_lights = AreaLightsImportanceSampler(self.area_lights)
        ggx = GGXSampler()
        self.direct_sampler = MISSampler(area_lights, ggx, 0.5)
        self.indirect_sampler = area_lights

        start = perf_counter()
        last_display = start
        self.russia_roulette_level = 1
        self.russian_roulette_prob = 1.0
        self.max_bounce = 1
        num = 50
        print (f"Rendering (  )...")
        print (f" - samples {num}")
        print (f" - bounces {self.max_bounce}")
        print (f" - bounce russian roulette starting from bounce {self.russia_roulette_level} with probability {self.russian_roulette_prob}")
        print (f" - using sampler {self.direct_sampler.get_name()} for direct lighting")
        print (f" - using sampler {self.indirect_sampler.get_name()} for indirect lighting")
        for y in range(h):
            for x in range (w):
                px = (2 * (x + 0.5) / w - 1) * htx
                py = (1 - 2 * (y + 0.5) / h) * hty
                d = (cam_forward + cam_right * px + cam_up * py).normalized()
                # if x != 130 or y != h-1-188: continue
                image[(h-1-y) * w + x] = self.primary(cam_pos, d, num)

                now = perf_counter()
                if now - last_display > 4.0:
                    last_display = now
                    ratio = (y * w + x) / (w*h)
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
            self.draw_data = MonteCarloIndirectDrawData(dimensions)

        self.draw_data.draw()

        self.unbind_display_space_shader()
        gpu.state.blend_set('NONE')


class MonteCarloIndirectDrawData:
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
    bpy.utils.register_class(MonterCarloIndirectEngine)

    for panel in get_panels():
        panel.COMPAT_ENGINES.add('MonteCarlo')


def unregister():
    bpy.utils.unregister_class(MonterCarloIndirectEngine)

    for panel in get_panels():
        if 'MonteCarlo' in panel.COMPAT_ENGINES:
            panel.COMPAT_ENGINES.remove('MonteCarlo')


if __name__ == "__main__":
    register()