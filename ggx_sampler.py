import bpy
from mathutils import Vector
import math
import random
from math import radians

# Adding primiteves

def add_cube(point,name, size=1.0):
    bpy.ops.mesh.primitive_cube_add(
        size=size,
        location=point if isinstance(point, tuple) else (point.x, point.y, point.z)
    )
    obj = bpy.context.active_object
    obj.name = name
    return obj

def add_sphere(point, name, mat = None,radius=1.0):
    bpy.ops.mesh.primitive_ico_sphere_add(
        radius=radius,
        location=point,
        subdivisions=6
    )
    if not mat is None:
        bpy.context.active_object.data.materials.append(mat)
    obj = bpy.context.active_object
    obj.name = name
    for poly in obj.data.polygons:
        poly.use_smooth = True
    
    return bpy.context.active_object

def add_text(point, name, text, scale = 0.1, rotation = Vector ((90, 0, -90)), mat = None):
    bpy.ops.object.text_add(location=point)
    text_obj = bpy.context.object
    text_obj.name = name
    text_obj.data.body =text
    text_obj.rotation_euler.x = radians(rotation.x)
    text_obj.rotation_euler.y = radians(rotation.y)
    text_obj.rotation_euler.z = radians(rotation.z)
    text_obj.scale = (scale, scale, scale)
    text_obj.data.align_x = 'CENTER'
    if not mat is None:
        text_obj.data.materials.append(mat)    
    return text_obj

def add_pdf_visualization(point, name, S):
    obj = add_sphere(point, name)
    mesh = obj.data
    for v in mesh.vertices:
        p = v.co.normalized()
        cds = math.sqrt(4*math.pi / len(mesh.vertices)) * 2.0
        value = 0
        for wi, pdf in S:
            if pdf == 0.0: continue
            ds = (wi-p).length_squared
            if ds < cds:
                value = pdf
                cds = ds
        v.co = p * value
    return obj

def add_f_visualization(point, name, S, f, div_pdf = False):
    obj = add_sphere(point, name)
    mesh = obj.data
    for v in mesh.vertices:
        p = v.co.normalized()
        cds = math.sqrt(4*math.pi / len(mesh.vertices)) * 2.0
        pdf_value = 0
        for wi, pdf in S:
            ds = (wi-p).length_squared
            if ds < cds:
                pdf_value = pdf
                cds = ds
        r = 0.0
        if pdf_value > 0.000001:
            r = f(p)
            if div_pdf: r /= pdf_value
        v.co = p * r
    return obj

# Project function on sphere

def hemi_project(f, obj):
    num_samples = 0
    mesh = obj.data
    for v in mesh.vertices:
        wi = v.co.normalized()
        r = 0.0
        if wi.z > 0.0:
            r = f(wi)
            num_samples += 1
        v.co = wi * r
    return num_samples

def sphere_project(f, obj):
    num_samples = 0
    mesh = obj.data
    for v in mesh.vertices:
        wi = v.co.normalized()
        r = f(wi)
        num_samples += 1
        v.co = wi * r
    return num_samples


# Samplers

def uniform_sample ():
    u = random.uniform(0.0, 2.0 * math.pi)
    cos_v = random.uniform(0.0, 1.0)
    v = math.acos(cos_v)

    x = math.sin(v) * math.cos(u)
    y = math.sin(v) * math.sin(u) 
    z = math.cos(v)

    pdf = 1.0 / (2.0 * math.pi)
    
    return Vector((x,y,z)), pdf

def pdf_uniform():
    return 1.0 / (2.0 * math.pi)

def cosine_sample ():
    u = random.uniform(0.0, 1.0)
    v = random.uniform(0.0, 2.0 * math.pi)

    r = math.sqrt(u)
    x = r*math.cos(v)
    y = r*math.sin(v)
    z = math.sqrt(1.0-u)

    pdf = z / math.pi

    return Vector((x,y,z)), pdf

def pdf_cosine(wi):
    return wi.z / math.pi

import math
from mathutils import Vector
import random

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
        return Vector((0, 0, 0)), 0.0

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
        return wi, 0.0

    # ---- PDF (VNDF-consistent) ----
    n_dot_h = max(0.0, n.dot(h))
    wo_dot_h = max(0.0, wo.dot(h))
    n_dot_wo = max(0.0, n.dot(wo))

    if wo_dot_h <= 0.0 or n_dot_wo <= 0.0:
        return wi, 0.0

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

    return wi, pdf


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

# BDRF definition

def lambert_BDRF(kd):
    return kd / math.pi

def phong_BRDF(wi, wr, N):
    k_s = 0.75
    shininess = 32.0
    
    R = 2.0 * N.dot(wi) * N - wi
    
    specular = k_s * math.pow(max(0, R.dot(wr)), shininess) #(shininess + 2)/(2*math.pi)
    
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


# Li

lights=[]
for obj in bpy.data.objects:
     if obj.name.startswith('LightProxy'):
         lights.append(obj)

def Li(wi, lights):
    s = Vector((0,0,0))
    e = wi*10.0
    for light in lights:
        ls = light.matrix_world.inverted() @ s
        le = light.matrix_world.inverted() @ e
        result, location_local, normal_local, face_index = light.ray_cast(ls,le)
        
        if result: return 1.0#max(0, (3.0 - (location_local-s).length)) / 3.0
    return 0.0 # ambient


# RIS

def RIS_sample(sampler, f, M):
    def select_weighted(W, total_weight):
        sw = random.uniform(0, total_weight)
        c = 0
        for i, w in enumerate(W):
            c += w
            if sw <= c: return i
        return len(W) - 1
    
    S = [sampler() for i in range(M)]
    W = [f(s) / pdf for (s, pdf) in S]
    total_weight = sum(W)

    i = select_weighted(W, total_weight)
    x = S[i][0]
    if W[i] == 0.0: return x, 0.0
    wx = (1.0/f(x)) * (total_weight / M)
    return x, wx


# Delete generated objects

bpy.ops.object.select_all(action='DESELECT')
for obj in bpy.data.objects:
     if obj.name.startswith('F_') or obj.name.startswith('S_'):
         obj.select_set(True)
bpy.ops.object.delete()


# Get the light and camera

cam = bpy.data.objects.get("Camera")
cam_matrix = cam.matrix_world
cam_pos = cam_matrix.to_translation()
wo = cam_pos.normalized()


# Monte Carlo
def monte_carlo_estimate(f, S):
    l = 0.0
    for wi, pdf in S:
        if pdf == 0.0: continue
        l += f(wi) / pdf
    return l / len(S)


# Display BDRF, Li

kd = 1.0
roughness=0.5
N = Vector((0,0,1))


# Display samplers
def div_pdf(f, pdf):
    if pdf < 0.0001: return 0.0
    return f / pdf
wo = Vector((-0.0069, -0.9858, -0.1677))
N = Vector((0.0000, -1.0000, 0.0000))
roughness = 0.800000011920929
sphere_project(lambda wi: ggx_pdf_vndf(wi, wo, N, roughness), add_sphere(Vector((0, 0, 0)), "F_GGX_PDF"))
sphere_project(lambda wi: ggx_BDRF(wi, wo, N, Vector((0.04, 0.04, 0.04)), roughness).length, add_sphere(Vector((0, 0, 0)), "F_GGX_BDRF"))
S = [ggx_sample_vndf(wo, N, roughness) for _ in range (15000)]
add_pdf_visualization(Vector((0, 0, 0)), "F_GGX_ACC", S)
print (ggx_BDRF(N, N, N, Vector((0.04, 0.04, 0.04)), roughness).length)
    # # cosine
    # sphere_project(lambda wi: div_pdf(l_f(wi), pdf_cosine(wi)), add_sphere(Vector((0, 0, 3)), "F_Cosine_Final"))
    # sphere_project(lambda wi: pdf_cosine(wi), add_sphere(Vector((0, 1, 3)), "F_Cosine_PDF"))
    # # uniform
    # sphere_project(lambda wi: l_f(wi) / pdf_uniform(), add_sphere(Vector((-4, 0, 3)), "F_Uniform_Final"))
    # sphere_project(lambda wi: pdf_uniform(), add_sphere(Vector((-4, 1, 3)), "F_Uniform_PDF"))



S = [ggx_sample_vndf(Vector((-0.0069, -0.9858, -0.1677)), Vector((0.0000, -1.0000, 0.0000)), 0.800000011920929) for i in range (30)]
for s in S:
    wi, pdf = s
    n_pdf = ggx_pdf_vndf(wi, Vector((-0.0069, -0.9858, -0.1677)), Vector((0.0000, -1.0000, 0.0000)), 0.800000011920929)
    print (pdf, n_pdf)
