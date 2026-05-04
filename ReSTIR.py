import sys
from enum import Enum
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

def luminance_rec2020(rgb):
    return 0.2627 * rgb.x + 0.6780 * rgb.y + 0.0593 * rgb.z

def get_cosmetic_duration(t):
    ms = int(t * 1000) % 1000
    t = int(t)

    d, t = divmod(t, 86400)
    h, t = divmod(t, 3600)
    m, s = divmod(t, 60)

    parts = [f"{v}{u}" for v, u in [(d,"d"), (h,"h"), (m,"m"), (s,"s")] if v]
    ret = " ".join(parts)
    if len(ret) == 0: return "0s"
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

def sample_disk(x, y, r):
    radius = r * math.sqrt(random.uniform(0, 1))
    angle = random.uniform(0, 2*math.pi)
    dx = int(x + math.cos(angle) * radius)
    dy = int(y + math.sin(angle) * radius)
    return dx, dy

# Area light =================

class AreaLight:
    obj = None
    matrix_world = None
    inv_matrix_world = None
    mat = None
    li = None
    verts = None
    ggx_mat_cache = None
    polygons = []
    area = None
    weight = None
    name = None

    def __init__(self, instance, ggx_mat_cache):
        self.obj = instance.object
        self.ggx_mat_cache = ggx_mat_cache
        obj = instance.object
        self.name = obj.name
        mesh = obj.data
        self.matrix_world = instance.matrix_world.copy()
        self.inv_matrix_world = instance.matrix_world.inverted().copy()
        self.verts = [instance.matrix_world @ v.co for v in mesh.vertices]
        self.mat = self.ggx_mat_cache.get(obj.material_slots[0].material)
        k_d, k_s, k_r, k_e, k_es = self.mat
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

    def set_params(self, param_target_pdf = None, param_mat = None):
        self.param_target_pdf = param_target_pdf
        self.param_mat = param_mat
        for s in self.samplers:
            s.set_params(param_target_pdf, param_mat)

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    @abstractmethod
    def samples(self, p, n, wo, num):
        pass
    
    # return (wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta)
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    def sample (self, p, n, wo):
        return self.samples(p, n, wo, 1)[0]
        
    # return pdf
    @abstractmethod
    def pdf(self, p, n, wo, wi):
        pass
    @abstractmethod
    def get_name(self):
        pass
    # checks sample for l_data
    def check_l_data(s, context_name):
        wi, l_data, pdf, mis_w, cos_theta = s
        if not l_data:
            print (f"{context_name} : Light data missing from sample {Sampler.debug(s)}")
            return False
        return True
    def get_pdf_rays(self):
        num = 0
        for s in self.samplers:
            num += s.get_pdf_rays ()
        return num
    # returns a dictionary with the content of the sample
    def debug(s):
        if s is None: return {}
        wi, l_data, pdf, mis_w, cos_theta = s
        l_data_dict = {}
        if l_data is not None:
            pos, nor, li, l = l_data
            l_data_dict = {"pos" : pos, "nor" : nor, "li": li, "l" : l.name if l is not None else "unknown"}
        return {"wi" : wi, "l_data" : l_data_dict, "pdf" : pdf, "mis_w" : mis_w, "cos_theta" : cos_theta}


class SimpleSampler(Sampler):
    sample_l = None
    pdf_l = None
    name = None

    def __init__(self, sample_l = cosine_sample, pdf_l = cosine_pdf, name = "Cosine"):
        self.sample_l = sample_l
        self.pdf_l = pdf_l
        self.name = name

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    def samples (self, p, n, wo, num):
        S = []
        for _ in range(num):
            wi, pdf, cos_theta = self.sample_l(n)
            S.append((wi, None, pdf, 1.0, cos_theta))
        return S
    
    def pdf(self, p, n, wo, wi):
        return self.pdf_l(p, n, wi)

    def get_name(self):
        return self.name

class AreaLightsImportanceSampler(Sampler):
    area_lights : List[AreaLight] = []
    weight = 0

    pdf_rays = 0

    def __init__(self, area_lights : List[AreaLight]):

        self.area_lights = area_lights
        for l in self.area_lights:
            self.weight += l.weight

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    def samples(self, p, n, wo, num):
        S = []
        for _ in range(num):
            l = self.area_lights[get_weighted_random_index(self.area_lights, self.weight, lambda l: l.weight)]
            pos, nor = l.get_weight_random_point()
            pdf, wi = l.pdf(p, n, pos, nor)
            pdf *= l.weight / self.weight # adjust pdf based on light selection pdf
            cos_theta = max(0, wi.dot(n))
            S.append((wi, (pos, nor, l.li, l), pdf, 1.0, cos_theta))
        return S
    
    def pdf(self, p, n, wo, wi):
        self.pdf_rays += 1
        pdf = 0.0
        for l in self.area_lights:
            hit, pos, nor = l.ray_trace(p, wi, SAMPLER_DISTANCE)
            if not hit: continue
            l_pdf, _ = l.pdf(p, n, pos, nor)
            pdf += l_pdf * l.weight / self.weight
        return pdf
    
    def get_pdf_rays(self):
        num = self.pdf_rays
        for s in self.samplers:
            num += s.get_pdf_rays ()
        return num

    
    def get_name(self):
        return f"LightsImportance({len(self.area_lights)} lights)"
    
class GGXSampler(Sampler):
    def __init__(self):
        super().__init__()

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    def samples(self, p, n, wo, num):
        S = []
        k_d, k_s, k_r, k_e, k_es = self.param_mat
        for _ in range(num):
            wi, pdf, cos_theta = ggx_sample_vndf(wo, n, k_r)
            S.append ((wi, None, pdf, 1.0, cos_theta))
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

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
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
            wi1, l_data1, pdf1, mis_w1, cos_theta1 = s
            pdf2 = s2.pdf(p, n, wo, wi1)
            w = 0 if pdf1 == 0.0 else (num1 * pdf1) / (num1 * pdf1 + num2 * pdf2)
            S.append((wi1, l_data1, pdf1, w, cos_theta1))
        for s in S2:
            wi2, l_data2, pdf2, mis_w2, cos_theta2 = s
            pdf1 = s1.pdf(p, n, wo, wi2)
            w = 0 if pdf2 == 0.0 else (num2 * pdf2) / (num1 * pdf1 + num2 * pdf2)
            S.append((wi2, l_data2, pdf2, w, cos_theta2))
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

    # return [(wi, l_data | None, pdf, mis_w(1.0 - if not MIS), cos_theta), ... ]
    # for area light importance samplers l (pos, nor, li, l) represents the sample pos, nor, li on light l - None otherwise
    def samples(self, p, n, wo, num):
        F = []

        s = self.samplers[0]
        for _ in range(num):
            S = s.samples(p, n, wo, self.M)
            W = [0.0 if pdf == 0.0 or mis_w == 0.0 else self.param_target_pdf(wi, cos_theta) * mis_w / pdf for wi, l_data, pdf, mis_w, cos_theta in S]
            total_weight = sum(W)

            champion = get_weighted_random_index(W, total_weight, lambda w: w)
            wi, l_data, pdf, mis_w, cos_theta = S[champion]
            if W[champion] == 0.0:
                F.append((wi, l_data, 0.0, 0.0, cos_theta))
            else:
                f = W[champion] * pdf / mis_w # compute target_pdf from weight instead of calculating it
                w = (total_weight / self.M) / f
                # RIS replaces the mis_w as 1 it's already incorporated into w. We can have RIS from MIS
                # and 1/w can be used as a PDF in the sense that f(x) / pdf = f(x) * w
                F.append((wi, l_data, 1/w, 1.0, cos_theta)) 

        return F
    
    def pdf(self, p, n, wo, wi):  # RIS Samples cannot be used in any strategy that requires computing the PDF for a generic wi
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
        return Vector((0, 0, 0))
    
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

def bdrf(wi, wo, n, m):
    k_d, k_s, k_r, k_e, k_es = m
    return lambert_BDRF(k_d) + ggx_BDRF (wi, wo, n, Vector((0.04, 0.04, 0.04)), k_r) * k_s

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

class Reservoir:
    s = None # (wi, l_data, pdf, mis_w, cos_theta)
    w_sum = 0
    c_sum = Vector((0, 0, 0))
    m = 0

    def __init__(self):
        self.s = None # (wi, l_data, pdf, mis_w, cos_theta)
        self.w_sum = 0
        self.c_sum = Vector((0, 0, 0))
        self.m = 0

    def read(self, t):
        self.s, self.w_sum, c_sum, self.m = t
        self.c_sum = c_sum.copy()
        return self
    
    def write (self):
        return (self.s, self.w_sum, self.c_sum.copy(), self.m)
    
    # s = sampler sample (wi, l_data, pdf, mis_w(1.0 - if not MIS), cos_theta)
    def add_sample(self, s, color):
        if not Sampler.check_l_data(s, "Reservoir.add_sample"): return

        w = luminance_rec2020(color)
        
        self.w_sum += w
        self.c_sum += color
        self.m += 1

        if random.random() * self.w_sum < w:
            self.s = s
        return self

    def add_reseroir(self, other):
        self.w_sum += other.w_sum
        self.c_sum += other.c_sum
        self.m += other.m
        if random.random() * self.w_sum <= other.w_sum:
            self.s = other.s
        return self

    def decay(self, sum_decay, m_decay):
        if self.m == 0: return self

        new_m = max (1, int(self.m * m_decay))
        m_factor = new_m / self.m
        self.m = new_m
        self.w_sum *= sum_decay * m_factor
        self.c_sum *= sum_decay * m_factor

        return self
    
    def recompute_pdf_cos_theta(self, p, n, wo, sampler : Sampler):
        if self.s is None: return self
        wi, l_data, pdf, mis_w, cos_theta = self.s
        if cos_theta == 0.0: return self
        pos, nor, li, l = l_data
        
        wi = (pos-p).normalized()
        new_cos_theta = max(0, n.dot(wi))
        new_pdf = sampler.pdf (p, n, wo, wi)
        # Some samplers such as MIS and RIS don't supoort recomputation of pdf
        new_pdf = pdf if new_pdf is None or new_pdf == 0.0 else new_pdf

        self.c_sum *= (new_cos_theta / cos_theta) * (pdf / new_pdf)
        self.w_sum = luminance_rec2020(self.c_sum)

        self.s = (wi, l_data, new_pdf, mis_w, new_cos_theta)
        
        return self

    def recompute_pdf_cos_theta_jacobian(self, p, n, d_p, wo, sampler : Sampler):
        if self.s is None: return self
        d_wi, l_data, pdf, mis_w, cos_theta = self.s
        pos, nor, li, l = l_data

        wi = (pos-p).normalized()
        new_cos_theta = max(0, n.dot(wi))

        # jacobian
        d_r_sq = (pos - d_p).length_squared
        d_cos_light = abs(nor.dot(-d_wi))

        r_sq = (pos - p).length_squared
        cos_light = abs(nor.dot(-wi))

        j = (cos_light / d_cos_light) * (d_r_sq / r_sq)
        
        new_pdf = sampler.pdf (p, n, wo, wi) * j

        self.c_sum *= (new_cos_theta / cos_theta) * (pdf / new_pdf)
        self.w_sum = luminance_rec2020(self.c_sum)
        self.s = (wi, l_data, new_pdf, mis_w, new_cos_theta)
        
        return self

    def debug(t):
        if t is None: return {}
        s, w_sum, c_sum, m = t
        s_dict = {}
        if s is not None:
            s_dict = Sampler.debug(s)
        return {"s" : s_dict,
                "w_sum" : w_sum,
                "c_sum" : c_sum,
                "m" : m}

class ReSTIREngine(bpy.types.RenderEngine):
    # These three members are used by Blender to set up the
    # RenderEngine; define its internal name, visible name and capabilities.
    bl_idname = "ReSTIR"
    bl_label = "ReSTIR"
    bl_use_preview = False
    bl_use_shading_nodes = True
    bl_use_world_space_shading = True

    area_lights : List[AreaLight] = []
    ggx_mat_cache : GGXMaterialCache = None
    sampler : Sampler = None
    gbuffer = []
    reservoirs = []

    # configs
    M = 10 # numbers of samples to be merged in reservoir initialization
    missing_reservoir_color = [0, 0, 0, 1]
    # when joining reservoirs (or shading) the reservoir sample pdf ^ cost_theta is computed and 
    # in combination with the orignal pdf & cos_theta the weight,color and sample are 
    # adjusted to reflect the new values. This makes the falloff of area lights to be properly
    # light rays are used to compute the pdf for area lights
    recompute_pdf_cos_theta = True
    spatial_radius = 5
    spatial_num = 5
    spatial_distance_threshold = 0.04
    spatial_nors_threshold = 0.9
    spatial_halfs_threshold = 0.83

    # stats
    rays = 0
    light_rays = 0

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
    
    def update_render_passes(self, scene=None, render_layer=None):
        print ("update_render_passes")
        self.register_pass(scene, render_layer, "init", 4, 'FLOAT', 'COLOR')
        self.register_pass(scene, render_layer, "shadow", 4, 'FLOAT', 'COLOR')
        self.register_pass(scene, render_layer, "temporal", 4, 'FLOAT', 'COLOR')
        self.register_pass(scene, render_layer, "spatial", 4, 'FLOAT', 'COLOR')

    
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
        self.rays += 1
        scene = bpy.context.scene
        depsgraph = bpy.context.evaluated_depsgraph_get()
        result, pos, nor, index, object, matrix = scene.ray_cast(depsgraph, s, wi, distance = d)
        if result:
            mesh = object.data
            mat = self.ggx_mat_cache.get(object.material_slots[mesh.polygons[index].material_index])
            return result, pos, nor, mat
        return False, None, None, None
    
    def visibility(self, s, e) -> float:
        wi = e-s
        d = wi.length
        wi /= d
        hit, pos, nor, mat = self.ray_trace(s + wi * SAMPLER_BIAS, wi, d)
        if not hit: return 1.0
        return 1.0 if (pos - e).length < 3.0 * SAMPLER_BIAS else 0.0
    
    def ray_trace_lights(self, p, wi, d):
        self.light_rays += 1
        ret_l = None
        ret_d = d
        ret_pos = None
        ret_nor = None
        for l in self.area_lights:
            hit, pos, nor = l.ray_trace(p, wi, ret_d)
            if not hit: continue
            ret_d = (pos-p).length
            ret_pos = pos
            ret_nor = nor
            ret_l = l
        return (True, (ret_pos, ret_nor, ret_l.li, ret_l)) if ret_l is not None else (False, None)
    
    def fill_samples_light_data (self, S, p, n):
        for i, s in enumerate(S):
            wi, l_data, pdf, mis_w, cos_theta = s
            if l_data is not None: continue
            hit, l_data = self.ray_trace_lights(p + n * SAMPLER_BIAS, wi, SAMPLER_DISTANCE)
            l_data = l_data if l_data is not None else (wi, Vector((0, 0, 0)), Vector ((0, 0, 0)), None)
            S[i] = (wi, l_data, pdf, mis_w, cos_theta)
   
    # Passes

    # params = cam_pos
    def generate_gbuffer(self, x, y, wo, params):
        cam_pos = params
        hit, pos, nor, mat = self.ray_trace(cam_pos, -wo, SAMPLER_DISTANCE)
        self.gbuffer[y*self.size_x+x] = (pos, nor, mat) if hit else None

    # params = destination reservoir
    def init_reservoir (self, x, y, wo, params):
        dst = params
        ofs = y*self.size_x+x

        # read gbuffer
        gbuffer = self.gbuffer[ofs]
        if not gbuffer: dst[ofs] = None; return
        p, n, m = gbuffer

        # skip emissive materials (area lights)
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0: return

        def f(s)->Vector:
            wi, l_data, pdf, mis_w, cos_theta = s
            pos, nor, li, l = l_data
            return bdrf(wi, wo, n, m) * li * cos_theta * mis_w / pdf
      
        r = Reservoir()
        self.sampler.set_params(param_mat=m, param_target_pdf=lambda wi, cos_theta: luminance_rec2020(bdrf(wi, wo, n, m)))
        S = self.sampler.samples(p, n, wo, self.M)
        self.fill_samples_light_data (S, p, n)
        for s in S:
            wi, l_data, pdf, mis_w, cos_theta = s
            if pdf == 0.0: continue
            r.add_sample(s, f(s))

        if r.s is None:dst[ofs] = None; return

        dst[ofs] = r.write()

    # params = reservoir to shadow
    def shadow(self, x, y, wo, params):
        dst = params
        ofs = y*self.size_x+x

        # read gbuffer
        gbuffer = self.gbuffer[ofs]
        if not gbuffer: dst[ofs] = None; return
        p, n, m = gbuffer

        # skip emissive materials (area lights)
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0: return

        r_data = dst[ofs]
        if r_data is None: return
        r = Reservoir().read(r_data)
        wi, l_data, pdf, mis_w, cos_theta = r.s
        pos, nor, li, l = l_data
        v = self.visibility(p, pos)
        r.w_sum *= v
        r.c_sum *= v
        dst[ofs] = r.write()

    def get_spatial_matching_reservoir(self, wo, p, n, h, m, dx, dy, src):
        if dx < 0 or dx >= self.size_x or dy < 0 or dy >= self.size_y: return None
        ofs = dy * self.size_x + dx

        # sample gbuffer
        gbuffer = self.gbuffer[ofs]
        if gbuffer is None: return None
        pos, nor, mat = gbuffer

        # avoid emissive (area lights)
        k_d, k_s, k_r, k_e, k_es = mat
        if k_es > 0: return None

        # sample reseroir
        r_data = src[ofs]
        if r_data is None: return None
        r = Reservoir().read(r_data)

        # Roughness similarity weight
        curr_k_d, curr_k_s, curr_k_r, curr_k_e, curr_k_es = m
        if max(0, 1.0 - abs(curr_k_r - k_r) * 2.0) < 0.2:
            return None
        # half vector
        r_wi, r_l_data, r_pdf, r_mis_w, r_cos_theta = r.s

        if h is not None:
            halfs_dot = h.dot((r_wi+wo).normalized()) - self.spatial_halfs_threshold
            if halfs_dot <= 0.0: return None

        # normal
        dot = n.dot(nor) - self.spatial_nors_threshold
        if dot <= 0: return None

        # position
        d = (p-pos).length
        if d > self.spatial_distance_threshold: return None

        return r

    # params = (src, dst)
    def spatial_reuse(self, x, y, wo, params):
        src, dst = params
        ofs = y*self.size_x+x

        # read gbuffer
        gbuffer = self.gbuffer[ofs]
        if not gbuffer: dst[ofs] = None; return
        p, n, m = gbuffer

        # skip emissive materials (area lights)
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0: dst[ofs] = None; return

        r = Reservoir()
        r_data = src[ofs]
        h = None
        if r_data is not None:
            r.read(r_data)
            # half vector
            wi, l_data, pdf, mis_w, cos_theta = r.s

            h = (wi + wo).normalized()

        for i in range (self.spatial_num):
            dx, dy = sample_disk(x, y, self.spatial_radius)
            dr = self.get_spatial_matching_reservoir(wo, p, n, h, m, dx, dy, src)
            if dr is None: continue
            if random.random() < self.spatial_shadowing_ratio:
                d_wi, d_l_data, d_pdf, d_mis_w, d_cos_theta = dr.s
                d_pos, d_nor, d_li, d_l = d_l_data
                v = self.visibility(p, d_pos)
                dr.w_sum *= v
                dr.c_sum *= v

            r.add_reseroir(dr.recompute_pdf_cos_theta(p, n, wo, self.sampler) if self.recompute_pdf_cos_theta else dr)

        dst[ofs] = None if r.s is None else r.write()

    # params (prev, dst)
    def temporal_reuse(self, x, y, wo, params):
        prev, dst = params
        ofs = y*self.size_x+x

        r_data = dst[ofs]
        if r_data is None: return

        prev_r_data = prev[ofs]
        if prev_r_data is None: return

        p_r = Reservoir().read(prev_r_data)
        r = Reservoir().read(r_data).add_reseroir(p_r.decay(0.8, 0.8))

        dst[ofs] = r.write()

    # params = (source reservoir, destination image)
    def shade(self, x, y, wo, params):
        src, dst = params

        src_ofs = y*self.size_x+x
        dst_ofst = (self.size_y - 1 - y) * self.size_x + x

        # read gbuffer
        gbuffer = self.gbuffer[src_ofs]
        if not gbuffer: dst[dst_ofst] = [0, 0, 0, 1]; return
        p, n, m = gbuffer

        # emissive materials (area lights)
        k_d, k_s, k_r, k_e, k_es = m
        if k_es != 0.0: 
            li = k_e * k_es
            dst[dst_ofst] = [li.x, li.y, li.z, 1];
            return
        
        r_data = src[src_ofs]
        if r_data is None: dst[dst_ofst] = self.missing_reservoir_color; return
        r = Reservoir().read(r_data)

        f = r.c_sum / r.m
        dst[dst_ofst] = [math.sqrt(f.x), math.sqrt(f.y), math.sqrt(f.z), 1]

    # This is the method called by Blender for both final renders (F12) and
    # small preview for materials, world and lights.
    def render(self, depsgraph):
        scene = depsgraph.scene
        scale = scene.render.resolution_percentage / 100.0
        w = self.size_x = int(scene.render.resolution_x * scale)
        h = self.size_y = int(scene.render.resolution_y * scale)

        color = [0.0, 0.0, 1.0, 1.0]
        image = [color] * self.size_x * self.size_y


        result = self.begin_result(0, 0, self.size_x, self.size_y)

        
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
        self.gbuffer = [None] * w * h
        self.reservoirs = [[None] * w * h for _ in range(3)] 
        random.seed(42)

        src = self.reservoirs[0]
        dst = self.reservoirs[1]
        prev = self.reservoirs[2]

        def do_pass (pass_function, params):
            start = perf_counter()
            last_display = start
            display = False
            check_count = 0
            for y in range(h):
                for x in range (w):
                    px = (2 * (x + 0.5) / w - 1) * htx
                    py = (1 - 2 * (y + 0.5) / h) * hty
                    wo = -(cam_forward + cam_right * px + cam_up * py).normalized()
                    pass_function(x, y, wo, params)

                    check_count += 1
                    if check_count > 40:
                        check_count = 0
                        now = perf_counter()
                        if display:
                            if now - last_display > 4.0:
                                last_display = now
                                r = (y*w+x) / (w*h)
                                print (f"\t\t{int(r * 100.0)}% in {get_cosmetic_duration(now - start)}")
                        else:
                            if now - start > 4.0:
                                display = True
                                print (f"")

            return perf_counter() - start

        self.sampler = MISSampler(AreaLightsImportanceSampler(self.area_lights), GGXSampler(), 0.5)
        T = 1
        self.M = 5
        self.missing_reservoir_color = [0, 0, 0, 1]#[0, 1, 1, 1.0]
        self.recompute_pdf_cos_theta = True
        # 5 <- spatial join percentage, specular & shadows improves but fireflies
        # 2.5 <- better but no - more fireflies
        self.spatial_radius = 10
        self.spatial_num = 12
        #0.04 <- 0.04 good shadows compared with Cycles
        #0.1 <- 0.04 acceptable shadows compared with Cycles
        #0.5 <- Needed to remove artefacts
        self.spatial_distance_threshold = 0.4
        self.spatial_nors_threshold = 0.9 # <- Needed to remove artefacts
        # 0.98 # <- Good specular compared with Cycles
        self.spatial_halfs_threshold = 0.83 # <- Good specular
        self.spatial_shadowing_ratio = 0.0

        print (f"Rendering (  )...")
        print (f"join rejection pos > {self.spatial_distance_threshold} nors > {self.spatial_nors_threshold} halfs > {self.spatial_halfs_threshold}")
        start = perf_counter()

        print (" - Generate GBuffer ...", end="", flush=True)
        duration = do_pass(self.generate_gbuffer, cam_pos)
        print (f"\t done in {get_cosmetic_duration(duration)}")
        self.rays = 0
        self.light_rays = 0

        def present(src, pass_name):
            do_pass(self.shade, (src, image))
            layer = result.layers[0].passes[pass_name]
            layer.rect = image
            self.update_result(result)

        for i in range(T):
            print (f" - {i+1}/{T} init {self.M} {self.sampler.get_name()} ...", end="", flush=True)
            duration = do_pass(self.init_reservoir, dst)
            print (f"\t done in {get_cosmetic_duration(duration)}")
            t = dst; dst = src; src = t;

            present(src, "init")

            print (f"\t- shadow ...", end="", flush=True)
            duration = do_pass(self.shadow, src)
            print (f"\t done in {get_cosmetic_duration(duration)}")

            present (src, "shadow")

            print (f"\t- temporal ...", end="", flush=True)
            duration = do_pass(self.temporal_reuse, (prev, src))
            print (f"\t done in {get_cosmetic_duration(duration)}")

            present (src, "temporal")

            print (f"\t- spatial radius {self.spatial_radius} samples {self.spatial_num} shadow ratio {self.spatial_shadowing_ratio} ...", end="", flush=True)
            duration = do_pass(self.spatial_reuse, (src, dst))
            print (f"\t done in {get_cosmetic_duration(duration)}")
            t = dst; dst = src; src = t;
            t = prev; prev = src; src = t;
        
            present(prev, "spatial")

        present (prev, "Combined")

        # stats
        print (f"- rays {int(self.rays / T)} SPP {self.rays / (T*w*h) : .2f}")
        print (f"- light rays {int(self.light_rays / T)} SPP {self.light_rays / (T*w*h) : .2f}")
        pdf_rays = self.sampler.get_pdf_rays()
        print (f"- pdf rays {int(pdf_rays / T)} SPP {pdf_rays / (T*w*h) : .2f}")
        confidence = 0
        num_confidence = 0
        for r_data in prev:
            if r_data is None: continue
            r = Reservoir().read(r_data)
            confidence += r.m
            num_confidence += 1
        confidence /= T * num_confidence
        print (f"Average M samples per pixel {confidence : .4f}")            

        print (f"Rendering done in {get_cosmetic_duration(perf_counter() - start)}")
        
        self.area_lights = []
        self.ggx_mat_cache = None

        # Here we write the pixel values to the RenderResult
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
            self.draw_data = ReSTIRDrawData(dimensions)

        self.draw_data.draw()

        self.unbind_display_space_shader()
        gpu.state.blend_set('NONE')


class ReSTIRDrawData:
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
    bpy.utils.register_class(ReSTIREngine)

    for panel in get_panels():
        panel.COMPAT_ENGINES.add('MonteCarlo')


def unregister():
    bpy.utils.unregister_class(ReSTIREngine)

    for panel in get_panels():
        if 'MonteCarlo' in panel.COMPAT_ENGINES:
            panel.COMPAT_ENGINES.remove('MonteCarlo')


if __name__ == "__main__":
    register()