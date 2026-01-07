#pragma once

#include <ap_fixed.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include <ap_int.h>

static const int N = 256;

// 16-bit fixed-point za re in im
typedef ap_fixed<16, 8, AP_RND_CONV, AP_SAT> fx_t;

struct cplx_t {
    fx_t re;
    fx_t im;
};

// AXI-Stream word: 32-bit payload (re[15:0], im[31:16])
// (no keep/strb/user/id/dest)
typedef ap_axis<32, 0, 0, 0> axis32_t;

// Top-level AXI-Stream FFT
void fft256(
    hls::stream<axis32_t> &s_in,
    hls::stream<axis32_t> &s_out
);
