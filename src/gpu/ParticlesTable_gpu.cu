#include "ParticlesTable.hpp"
#include <cuda_runtime.h> 
#include <cstdio>

#define H2D cudaMemcpyHostToDevice
#define D2H cudaMemcpyDeviceToHost

#define CUDA_CHECK(call)                                                     \
    do {                                                                     \
        cudaError_t err = call;                                              \
        if (err != cudaSuccess) {                                            \
            fprintf(stderr, "CUDA error %s (%d) at %s:%d\n",                 \
                    cudaGetErrorString(err), err, __FILE__, __LINE__);       \
            std::abort();                                                    \
        }                                                                    \
    } while (0)

// Device-only function
// --- Direct N-body  ---
__global__ void dirnbody_kernel(const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict z,
                                const float* __restrict h,
                                const float* __restrict m,
                                float* __restrict ax,
                                float* __restrict ay,
                                float* __restrict az,
                                float* __restrict U,
                                int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    float xi = x[i], yi = y[i], zi = z[i], hi = h[i];
    float axt = 0.f, ayt = 0.f, azt = 0.f, Ut = 0.f;

    for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        float dx = x[j] - xi;
        float dy = y[j] - yi;
        float dz = z[j] - zi;
        float dr2 = dx*dx + dy*dy + dz*dz + 0.5f*(hi*hi + h[j]*h[j]);
        float invr = rsqrtf(dr2);       
        float invr3 = invr * invr * invr;

        float mjinvr3 = m[j] * invr3;
        axt += mjinvr3 * dx;
        ayt += mjinvr3 * dy;
        azt += mjinvr3 * dz;
        Ut  -= m[i] * m[j] * invr;
    }
    ax[i] = axt;  ay[i] = ayt;  az[i] = azt;  U[i] = Ut;
}

__global__ void dirnbody_2D_kernel(const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict h,
                                const float* __restrict m,
                                float* __restrict ax,
                                float* __restrict ay,
                                float* __restrict U,
                                int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    float xi = x[i], yi = y[i], hi = h[i];
    float axt = 0.f, ayt = 0.f, Ut = 0.f;

    for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        float dx = x[j] - xi;
        float dy = y[j] - yi;
        float dr2 = dx*dx + dy*dy + 0.5f*(hi*hi + h[j]*h[j]);
        float invr = rsqrtf(dr2);       
        float invr3 = invr * invr * invr;

        float mjinvr3 = m[j] * invr3;
        axt += mjinvr3 * dx;
        ayt += mjinvr3 * dy;
        Ut  -= m[i] * m[j] * invr;
    }
    ax[i] = axt;  ay[i] = ayt;  U[i] = Ut;
}

// --- BH Tree ---
__global__ void BHtree_kernel(const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict z,
                                const float* __restrict h,
                                const float* __restrict m,
                                float* __restrict ax,
                                float* __restrict ay,
                                float* __restrict az,
                                float* __restrict U,
                                const OctNode* __restrict nodes,
                                const int* __restrict order,
                                int root_idx,
                                float theta,
                                int N)
{
    // Get idx
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // Construct a stack to mimic the recursion
    const int MAX_STACK = 256;
    int stack[MAX_STACK];
    int sp = 0;                     // Current number of stack
    stack[sp++] = root_idx;

    // Get current quantities
    float xi = x[idx];
    float yi = y[idx];
    float zi = z[idx];
    float hi = h[idx];
    float mi = m[idx];

    float axi = 0.0f, ayi = 0.0f, azi = 0.0f, Ui = 0.0f;

    while (sp > 0) {
        int nidx = stack[--sp];  // Get node index
        const OctNode& node = nodes[nidx];

        if (node.pcount == 0) continue;

        if (node.isLeaf()) {
            for (int p = 0; p < node.pcount; ++p) {
                int j = order[node.ParticlesLocateidx + p];
                if (j == idx) continue;

                float dx = x[j] - xi;
                float dy = y[j] - yi;
                float dz = z[j] - zi;
                float r2 = dx*dx + dy*dy + dz*dz + 0.5f * (hi*hi + h[j]*h[j]);
                float invr = rsqrtf(r2);
                float invr3 = invr * invr * invr;
                axi += m[j] * dx * invr3;
                ayi += m[j] * dy * invr3;
                azi += m[j] * dz * invr3;
                Ui -= mi * m[j] * invr;
            }
            continue;
        }

        float dx = node.COMx - xi;
        float dy = node.COMy - yi;
        float dz = node.COMz - zi;
        float r2 = dx*dx + dy*dy + dz*dz;

        if (r2 > 1e-8f && node.cellsize() * rsqrtf(r2) < theta) {
            r2 += hi*hi;
            float invr = rsqrtf(r2);
            float invr3 = invr * invr * invr;
            float mInvr3 = node.Mtot * invr3;
            axi += mInvr3 * dx;
            ayi += mInvr3 * dy;
            azi += mInvr3 * dz;
            Ui -= mi * node.Mtot * invr;
        } else {
            for (int q = 0; q < 4; ++q) {
                int cidx = node.children[q];
                if (cidx >= 0 && sp < MAX_STACK){
                    stack[sp++] = cidx;
                }
            }
        }
    }
    // Store the acc into GPU array
    ax[idx] = axi;
    ay[idx] = ayi;
    az[idx] = azi;
    U[idx]  = Ui;
}

__global__ void BHtree_2D_kernel(const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict h,
                                const float* __restrict m,
                                float* __restrict ax,
                                float* __restrict ay,
                                float* __restrict U,
                                const QuadNode* __restrict nodes,
                                const int* __restrict order,
                                int root_idx,
                                float theta,
                                int N)
{
    // Get idx
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    // Construct a stack to mimic the recursion
    const int MAX_STACK = 128;
    int stack[MAX_STACK];
    int sp = 0;                     // Current number of stack
    stack[sp++] = root_idx;

    // Get current quantities
    float xi = x[idx];
    float yi = y[idx];
    float hi = h[idx];
    float mi = m[idx];

    float axi = 0.0f, ayi = 0.0f, Ui = 0.0f;

    while (sp > 0) {
        int nidx = stack[--sp];  // Get node index
        const QuadNode& node = nodes[nidx];

        if (node.pcount == 0) continue;

        if (node.isLeaf()) {
            for (int p = 0; p < node.pcount; ++p) {
                int j = order[node.ParticlesLocateidx + p];
                if (j == idx) continue;

                float dx = x[j] - xi;
                float dy = y[j] - yi;
                float r2 = dx*dx + dy*dy + 0.5f * (hi*hi + h[j]*h[j]);
                float invr = rsqrtf(r2);
                float invr3 = invr * invr * invr;
                axi += m[j] * dx * invr3;
                ayi += m[j] * dy * invr3;
                Ui -= mi * m[j] * invr;
            }
            continue;
        }

        float dx = node.COMx - xi;
        float dy = node.COMy - yi;
        float r2 = dx*dx + dy*dy;

        if (r2 > 1e-8f && node.cellsize() * rsqrtf(r2) < theta) {
            r2 += hi*hi;
            float invr = rsqrtf(r2);
            float invr3 = invr * invr * invr;
            float mInvr3 = node.Mtot * invr3;
            axi += mInvr3 * dx;
            ayi += mInvr3 * dy;
            Ui -= mi * node.Mtot * invr;
        } else {
            for (int q = 0; q < 4; ++q) {
                int cidx = node.children[q];
                if (cidx >= 0 && sp < MAX_STACK){
                    stack[sp++] = cidx;
                }
            }
        }
    }
    // Store the acc into GPU array
    ax[idx] = axi;
    ay[idx] = ayi;
    U[idx]  = Ui;
}



// --- Kick ---
__global__ void kick_kernel(float* __restrict vx,
                            float* __restrict vy,
                            float* __restrict vz,
                            const float* __restrict ax,
                            const float* __restrict ay,
                            const float* __restrict az,
                            const float* __restrict dt,
                            float scale,
                            int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float dt_i = dt[i];
    float fac  = scale * dt_i;
    vx[i] += fac * ax[i];
    vy[i] += fac * ay[i];
    vz[i] += fac * az[i];
}

__global__ void kick_2D_kernel( float* __restrict vx,
                                float* __restrict vy,
                                const float* __restrict ax,
                                const float* __restrict ay,
                                const float* __restrict dt,
                                float scale,
                                int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float dt_i = dt[i];
    float fac  = scale * dt_i;
    vx[i] += fac * ax[i];
    vy[i] += fac * ay[i];
}

// --- Drift ---
__global__ void drift_kernel(float* __restrict x,
                             float* __restrict y,
                             float* __restrict z,
                             const float* __restrict vx,
                             const float* __restrict vy,
                             const float* __restrict vz,
                             const float* __restrict dt,
                             float scale,
                             int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float dt_i = dt[i];
    float fac  = scale * dt_i;
    x[i] += fac * vx[i];
    y[i] += fac * vy[i];
    z[i] += fac * vz[i];
}

__global__ void drift_2D_kernel(float* __restrict x,
                                float* __restrict y,
                                const float* __restrict vx,
                                const float* __restrict vy,
                                const float* __restrict dt,
                                float scale,
                                int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    float dt_i = dt[i];
    float fac  = scale * dt_i;
    x[i] += fac * vx[i];
    y[i] += fac * vy[i];
}




// Method Definition
void ParticlesTable::device_init() {
    if (gpu_init) return;
    size_t n = static_cast<size_t>(N) * sizeof(float);

    CUDA_CHECK(cudaMalloc(&d_x , n)); CUDA_CHECK(cudaMalloc(&d_y , n)); CUDA_CHECK(cudaMalloc(&d_z , n));
    CUDA_CHECK(cudaMalloc(&d_vx, n)); CUDA_CHECK(cudaMalloc(&d_vy, n)); CUDA_CHECK(cudaMalloc(&d_vz, n));
    CUDA_CHECK(cudaMalloc(&d_m , n)); CUDA_CHECK(cudaMalloc(&d_h , n)); CUDA_CHECK(cudaMalloc(&d_dt, n));
    CUDA_CHECK(cudaMalloc(&d_ax, n)); CUDA_CHECK(cudaMalloc(&d_ay, n)); CUDA_CHECK(cudaMalloc(&d_az, n));
    CUDA_CHECK(cudaMalloc(&d_U , n));

    upload_all();  
    gpu_init = true;
}

void ParticlesTable::device_finalize() {
    if (!gpu_init) return;
    cudaFree(d_x);  cudaFree(d_y);  cudaFree(d_z);
    cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
    cudaFree(d_m);  cudaFree(d_h); cudaFree(d_dt);
    cudaFree(d_ax); cudaFree(d_ay); cudaFree(d_az); cudaFree(d_U);
    cudaFree(d_nodes_2D); cudaFree(d_nodes_3D); cudaFree(d_order);
    gpu_init = false;
}

void ParticlesTable::upload_all() {
    size_t n=N*sizeof(float);
    cudaMemcpyAsync(d_x , x.data() , n, H2D, stream);
    cudaMemcpyAsync(d_y , y.data() , n, H2D, stream);
    cudaMemcpyAsync(d_z , z.data() , n, H2D, stream);
    cudaMemcpyAsync(d_vx, vx.data(), n, H2D, stream);
    cudaMemcpyAsync(d_vy, vy.data(), n, H2D, stream);
    cudaMemcpyAsync(d_vz, vz.data(), n, H2D, stream);
    cudaMemcpyAsync(d_m , m.data() , n, H2D, stream);
    cudaMemcpyAsync(d_h , h.data() , n, H2D, stream);
    cudaMemcpyAsync(d_dt, dt.data(), n, H2D, stream);
    cudaStreamSynchronize(stream);
}

void ParticlesTable::download_state() {           
    size_t n=N*sizeof(float);
    cudaMemcpyAsync(x.data() , d_x , n, D2H, stream);
    cudaMemcpyAsync(y.data() , d_y , n, D2H, stream);
    cudaMemcpyAsync(z.data() , d_z , n, D2H, stream);
    cudaMemcpyAsync(vx.data(), d_vx, n, D2H, stream);
    cudaMemcpyAsync(vy.data(), d_vy, n, D2H, stream);
    cudaMemcpyAsync(vz.data(), d_vz, n, D2H, stream);
    cudaMemcpyAsync(_ax.data(),d_ax, n, D2H, stream);
    cudaMemcpyAsync(_ay.data(),d_ay, n, D2H, stream);
    cudaMemcpyAsync(_az.data(),d_az, n, D2H, stream);
    cudaMemcpyAsync(_U.data(), d_U,  n, D2H, stream);

    cudaStreamSynchronize(stream);
}

void ParticlesTable::calculate_a_dirnbody_gpu() {
    if (!gpu_init) device_init();
    int grid = (N + block - 1) / block;
    dirnbody_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y, d_z, d_h, d_m,
        d_ax, d_ay, d_az, d_U, N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::calculate_a_dirnbody_2D_gpu() {
    if (!gpu_init) device_init();
    int grid = (N + block - 1) / block;
    dirnbody_2D_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y, d_h, d_m,
        d_ax, d_ay, d_U, N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::calculate_a_BHtree_gpu() {
    // Prevent Tree recorvering
    cudaStreamSynchronize(stream);
    
    // Build quadtree
    OctTree tree = buildOctTree();

    // Get root
    const int root = tree.root_idx;

    // Get theta
    const float theta = tree.bhtheta;

    // Preparing sending tree to GPU
    int n_nodes = tree.nodes_list.size();
    int n_order = tree.order.size();

    // Allocate GPU
    // Free pointer first
    if (d_nodes_3D) cudaFree(d_nodes_3D);
    if (d_order)    cudaFree(d_order);  
    // Nodes
    CUDA_CHECK(cudaMalloc(&d_nodes_3D , n_nodes * sizeof(OctNode)));
    // Order
    CUDA_CHECK(cudaMalloc(&d_order ,  n_order * sizeof(int)));

    // Upload tree
    cudaMemcpyAsync(d_nodes_3D , tree.nodes_list.data() ,  n_nodes * sizeof(OctNode), H2D, stream);
    cudaMemcpyAsync(d_order , tree.order.data() ,  n_order * sizeof(int), H2D, stream);



    if (!gpu_init) device_init();
    int grid = (N + block - 1) / block;

    BHtree_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y,  d_z, d_h, d_m,
        d_ax, d_ay, d_az, d_U, d_nodes_3D, d_order,root, theta,  N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::calculate_a_BHtree_2D_gpu() {
    // Prevent Tree recorvering
    cudaStreamSynchronize(stream);

    // Build quadtree
    QuadTree tree = buildQuadTree();

    // Get root
    const int root = tree.root_idx;

    // Get theta
    const float theta = tree.bhtheta;

    // Preparing sending tree to GPU
    int n_nodes = tree.nodes_list.size();
    int n_order = tree.order.size();

    // Allocate GPU
    // Free pointer first
    if (d_nodes_2D) cudaFree(d_nodes_2D);
    if (d_order)    cudaFree(d_order);  
    // Nodes
    CUDA_CHECK(cudaMalloc(&d_nodes_2D , n_nodes * sizeof(QuadNode)));
    // Order
    CUDA_CHECK(cudaMalloc(&d_order ,  n_order * sizeof(int)));

    // Upload tree
    cudaMemcpyAsync(d_nodes_2D , tree.nodes_list.data() ,  n_nodes * sizeof(QuadNode), H2D, stream);
    cudaMemcpyAsync(d_order , tree.order.data() ,  n_order * sizeof(int), H2D, stream);



    if (!gpu_init) device_init();
    int grid = (N + block - 1) / block;

    BHtree_2D_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y, d_h, d_m,
        d_ax, d_ay, d_U, d_nodes_2D, d_order,root, theta,  N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::kick_gpu(float scale){
    int grid = (N + block - 1) / block;
    kick_kernel<<<grid, block, 0, stream>>>(
        d_vx, d_vy, d_vz,
        d_ax, d_ay, d_az,
        d_dt, scale, N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::drift_gpu(float scale){
    int grid = (N + block - 1) / block;
    drift_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y, d_z,
        d_vx, d_vy, d_vz,
        d_dt, scale, N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::kick_2D_gpu(float scale){
    int grid = (N + block - 1) / block;
    kick_2D_kernel<<<grid, block, 0, stream>>>(
        d_vx, d_vy,
        d_ax, d_ay,
        d_dt, scale, N);
    CUDA_CHECK(cudaGetLastError());
}

void ParticlesTable::drift_2D_gpu(float scale){
    int grid = (N + block - 1) / block;
    drift_2D_kernel<<<grid, block, 0, stream>>>(
        d_x, d_y,
        d_vx, d_vy,
        d_dt, scale, N);
    CUDA_CHECK(cudaGetLastError());
}