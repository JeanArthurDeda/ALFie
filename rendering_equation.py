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

def ggx_sample(wo, roughness):
    # Generate random numbers
    u1 = random.random()
    u2 = random.random()
    
    # Get alpha (roughness parameter)
    alpha = roughness * roughness
    
    # Sample microfacet normal (half-vector) using GGX distribution
    # This is the VNDF (Visible Normal Distribution Function) method
    # which samples visible microfacets for better efficiency
    
    # Step 1: Sample a random direction in the hemisphere
    phi = 2.0 * math.pi * u1
    
    # Sample slope distribution (Heitz 2018 method)
    # For GGX, we sample slopes directly
    a = 1.0 / alpha
    
    # Generate random point on disk with GGX distribution
    r = math.sqrt(u2 / (1.0 - u2))  # This gives distribution proportional to D()
    r = r * alpha  # Scale by roughness
    
    # Convert to slope (xy components)
    slope_x = r * math.cos(phi)
    slope_y = r * math.sin(phi)
    
    # Construct microfacet normal m from slopes
    # m = (-slope_x, -slope_y, 1).normalized()
    # But careful: this gives m in hemisphere around (0,0,1)
    inv_len = 1.0 / math.sqrt(slope_x * slope_x + slope_y * slope_y + 1.0)
    m = Vector((
        -slope_x * inv_len,
        -slope_y * inv_len,
        inv_len
    ))
    
    # Step 2: Reflect wo across m to get wi
    # wi = 2 * (wo · m) * m - wo
    wo_dot_m = wo.dot(m)
    wi = 2.0 * wo_dot_m * m - wo
    
    # Normalize (should be unit length already, but for safety)
    wi.normalize()
    
    # Step 3: Compute PDF
    pdf = pdf_ggx(wo, wi, roughness)
    
    return wi, pdf


def pdf_ggx(wi, wo, roughness):
    # Ensure directions are normalized
    wo = wo.normalized()
    wi = wi.normalized()
    
    # Normal is (0, 0, 1)
    N = Vector((0, 0, 1))
    
    # Compute half-vector
    # h = (wo + wi).normalized()
    h = wo + wi
    if h.length > 1e-6:
        h.normalize()
    else:
        return 0.0
    
    # Get alpha
    alpha = roughness * roughness
    alpha2 = alpha * alpha
    
    # Compute dot products
    NdotH = max(0.001, N.dot(h))
    NdotWo = max(0.001, N.dot(wo))
    NdotWi = max(0.001, N.dot(wi))
    HdotWi = max(0.001, h.dot(wi))
    
    # GGX NDF (D term)
    temp = (NdotH * NdotH) * (alpha2 - 1.0) + 1.0
    D = alpha2 / (math.pi * temp * temp)
    
    # Jacobian for half-vector to direction transformation
    # |dwi / dh| = 1 / (4 * (wi · h))
    # But for VNDF sampling, the PDF is:
    # pdf = D * NdotH / (4 * HdotWi)
    # This is the PDF of sampling wi using visible normal sampling
    
    # PDF of sampling half-vector
    pdf_h = D * NdotH
    
    # Transform to PDF of wi
    pdf = pdf_h / (4.0 * HdotWi)
    
    # Clamp to reasonable range
    pdf = max(0.0, min(1e6, pdf))
    
    return pdf


# BDRF definition

def lambert_BDRF(kd):
    return kd / math.pi

def phong_BRDF(wi, wr, N):
    k_d = 1.0
    diffuse = k_d * max(0, wi.dot(N)) / math.pi
    
    k_s = 0.75
    shininess = 32.0
    
    R = 2.0 * N.dot(wi) * N - wi
    
    specular = k_s * math.pow(max(0, R.dot(wr)), shininess) #(shininess + 2)/(2*math.pi)
    
    return diffuse + specular

def ggx_BDRF(wi, wo, N, F0, roughness):
    # Perfect mirror case (roughness = 0)
    if roughness < 1e-6:  # Treat as effectively zero
        # Compute perfect reflection direction
        R = 2.0 * N.dot(wi) * N - wi
        
        # Check if view direction aligns with reflection direction
        if abs(R.dot(wo) - 1.0) < 1e-6:
            # Perfect reflection - return Dirac delta represented as a large value
            # The actual BRDF is infinite, but in practice we return a very high value
            # and let the renderer handle it through importance sampling
            return 1e6  # Effectively infinite for numerical purposes
        else:
            return 0.0  # No reflection in other directions
    
    # Half-vector
    H = (wi + wo)
    H = H / H.length
    
    # Dot products with safety margin
    NdotWi = max(0.001, N.dot(wi))
    NdotWr = max(0.001, N.dot(wo))
    NdotH = max(0.001, N.dot(H))
    HdotWi = max(0.001, H.dot(wi))
    
    # 1. Normal Distribution Function (GGX)
    alpha = roughness * roughness
    alpha2 = alpha * alpha
    
    # Safe denominator for NDF
    temp = (NdotH * NdotH) * (alpha2 - 1.0) + 1.0
    denominator_D = math.pi * temp * temp
    D = alpha2 / denominator_D
    
    # 2. Geometry Function (Smith-GGX)
    def G1(NdotV):
        k = alpha / 2.0  # For Smith-GGX (k = α/2 for direct lighting)
        return NdotV / (NdotV * (1.0 - k) + k)
    
    G = G1(NdotWi) * G1(NdotWr)
    
    # 3. Fresnel (Schlick approximation)
    F = F0 + (1.0 - F0) * math.pow(1.0 - HdotWi, 5)
    
    # 4. Complete microfacet BRDF
    denominator = 4.0 * NdotWi * NdotWr
    specular = (F * D * G) / denominator
    
    return specular


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


# Display BDRF, Li and Rendering equation

ka = 1.0
roughness=0.5
N = Vector((0,0,1))

l_bdrf = lambda wi: lambert_BDRF(ka) + ggx_BDRF(wi, wo, N, 0.03, roughness=roughness)
l_li = lambda wi: Li(wi, lights)
l_f = lambda wi: l_bdrf(wi)*l_li(wi)*max(0, wi.z) # this is also the rendering equation
hemi_project(l_bdrf, add_sphere(Vector((0, 0, 0)), "F_BDRF"))
hemi_project(l_li, add_sphere(Vector((0, 0, 0)), "F_Li", mat=lights[0].material_slots[0].material))
hemi_project(l_f, add_sphere(Vector((0, 1.0, 0)), "F_F", mat= bpy.data.objects.get("mat_pdf").material_slots[0].material))


# Display samplers
# GGX
def div_pdf(f, pdf):
    if pdf < 0.0001: return 0.0
    return f / pdf
if False:
    sphere_project(lambda wi: div_pdf(l_f(wi), pdf_ggx(wi, wo, roughness)), add_sphere(Vector((4, 0, 3)), "F_GGX_Final"))
    sphere_project(lambda wi: pdf_ggx(wi, wo, roughness), add_sphere(Vector((4, 1, 3)), "F_GGX_PDF"))
    # cosine
    sphere_project(lambda wi: div_pdf(l_f(wi), pdf_cosine(wi)), add_sphere(Vector((0, 0, 3)), "F_Cosine_Final"))
    sphere_project(lambda wi: pdf_cosine(wi), add_sphere(Vector((0, 1, 3)), "F_Cosine_PDF"))
    # uniform
    sphere_project(lambda wi: l_f(wi) / pdf_uniform(), add_sphere(Vector((-4, 0, 3)), "F_Uniform_Final"))
    sphere_project(lambda wi: pdf_uniform(), add_sphere(Vector((-4, 1, 3)), "F_Uniform_PDF"))




num_samples = 30
num_ris_samples = 3
U = [uniform_sample() for _ in range(num_samples)]
C = [cosine_sample() for _ in range(num_samples)]
GGX = [ggx_sample(wo, roughness) for _ in range(num_samples)]
RIS = [(wi, 1/w if w != 0 else 0) for wi, w in [RIS_sample(cosine_sample, l_f, num_ris_samples) for _ in range(num_samples)]]
add_f_visualization(Vector((0, 1.5, 0)), "F_F_U", U, l_f)
add_f_visualization(Vector((0, 2, 0)), "F_F_C", C, l_f)
add_f_visualization(Vector((0, 2.5, 0)), "F_F_GGX", GGX, l_f)
add_f_visualization(Vector((0, 3, 0)), "F_F_RIS", RIS, l_f)


gt_value = monte_carlo_estimate(l_f, [cosine_sample() for i in range (16000)])
average_count = 50
u_value = 0.0
c_value = 0.0
ggx_value = 0.0
ris_value = 0.0
print ("==================================")
for i in range (average_count):
    U = [uniform_sample() for _ in range(num_samples)]
    C = [cosine_sample() for _ in range(num_samples)]
    GGX = [ggx_sample(wo, roughness) for _ in range(num_samples)]
    RIS = [(wi, 1/w if w != 0 else 0) for wi, w in [RIS_sample(cosine_sample, l_f, num_ris_samples) for _ in range(num_samples)]]
    u_value += monte_carlo_estimate(l_f, U)
    c_value += monte_carlo_estimate(l_f, C)
    ggx_value += monte_carlo_estimate(l_f, GGX)
    ris_value += monte_carlo_estimate(l_f, RIS)

u_value /= average_count
c_value /= average_count
ggx_value /= average_count
ris_value /= average_count

print (f"\Final")
print (f"u_value: {u_value}")
print (f"c_value: {c_value}")
print (f"ggx_value: {ggx_value}")
print (f"ris_value: {ris_value}")

add_text(Vector((-0.10, -0.2 + 1, 0.35)), "F_R_GT", f"GT\n{int(gt_value*10000)}", mat = lights[0].material_slots[0].material)
add_text(Vector((-0.10, -0.2 + 1.5,0.35)), "F_R_U", f"U\n{int(u_value*10000)}", mat = lights[0].material_slots[0].material)
add_text(Vector((-0.10, -0.2 + 2.0,0.35)), "F_R_C", f"COS\n{int(c_value*10000)}", mat = lights[0].material_slots[0].material)
add_text(Vector((-0.10, -0.2 + 2.5,0.35)), "F_R_GGX", f"GGX\n{int(ggx_value*10000)}", mat = lights[0].material_slots[0].material)
add_text(Vector((-0.10, -0.2 + 3.0,0.35)), "F_R_RIS", f"RIS({num_ris_samples})\n{int(ris_value*10000)}", mat = lights[0].material_slots[0].material)
add_text(Vector((-0.10, 2, 0.5)), "F_Desc", f"Aproximation of rendering equation with {num_samples} samples", mat = lights[0].material_slots[0].material)

