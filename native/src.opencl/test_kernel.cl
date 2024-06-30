__kernel void vec_add(__global const float *v1, __global const float *v2,
                             __local float *shared_mem, __global float *result)
{
    int wgid = get_group_id(0);
    int wgsize = get_local_size(0);
    int lid = get_local_id(0);

    shared_mem[lid] = v1[wgid * wgsize + lid] + v2[wgid * wgsize + lid];

    work_group_barrier(CLK_LOCAL_MEM_FENCE);

    int gid = get_global_id(0);
    result[gid] = shared_mem[lid];
}