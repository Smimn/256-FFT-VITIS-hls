#include "fft.h"

/*
 * Quarter-wave cosine tabela (0..90°)
 * ----------------------------------
 * Vsebuje cos(θ) za θ od 0° do 90° (65 točk).
 * FFT potrebuje sin/cos (twiddle faktorje).
 * Namesto računanja cos/sin v runtime uporabimo LUT
 * in simetrije kroga za generiranje vseh kotov.
 */
static const fx_t COS_Q[65] = {
    (fx_t)1.0000000000, (fx_t)0.9996988187, (fx_t)0.9987954562, (fx_t)0.9972904567,
    (fx_t)0.9951847267, (fx_t)0.9924795346, (fx_t)0.9891765100, (fx_t)0.9852776424,
    (fx_t)0.9807852804, (fx_t)0.9757021300, (fx_t)0.9700312532, (fx_t)0.9637760658,
    (fx_t)0.9569403357, (fx_t)0.9495281806, (fx_t)0.9415440652, (fx_t)0.9329927988,
    (fx_t)0.9238795325, (fx_t)0.9142097557, (fx_t)0.9039892931, (fx_t)0.8932243012,
    (fx_t)0.8819212643, (fx_t)0.8700869911, (fx_t)0.8577286100, (fx_t)0.8448535652,
    (fx_t)0.8314696123, (fx_t)0.8175848132, (fx_t)0.8032075315, (fx_t)0.7883464276,
    (fx_t)0.7730104534, (fx_t)0.7572088465, (fx_t)0.7409511254, (fx_t)0.7242470830,
    (fx_t)0.7071067812, (fx_t)0.6895405447, (fx_t)0.6715589548, (fx_t)0.6531728430,
    (fx_t)0.6343932842, (fx_t)0.6152315906, (fx_t)0.5956993045, (fx_t)0.5758081914,
    (fx_t)0.5555702330, (fx_t)0.5349976199, (fx_t)0.5141027442, (fx_t)0.4928981922,
    (fx_t)0.4713967368, (fx_t)0.4496113297, (fx_t)0.4275550934, (fx_t)0.4052413140,
    (fx_t)0.3826834324, (fx_t)0.3598950365, (fx_t)0.3368898534, (fx_t)0.3136817404,
    (fx_t)0.2902846773, (fx_t)0.2667127575, (fx_t)0.2429801799, (fx_t)0.2191012402,
    (fx_t)0.1950903220, (fx_t)0.1709618888, (fx_t)0.1467304745, (fx_t)0.1224106752,
    (fx_t)0.0980171403, (fx_t)0.0735645636, (fx_t)0.0490676743, (fx_t)0.0245412285,
    (fx_t)0.0000000000
};

/*
 * bitrev8
 * -------
 * Izvede bit-reversal za 8-bitno število (0..255).
 * FFT radix-2 zahteva, da so vhodni vzorci
 * najprej razporejeni v bit-reversed vrstni red.
 *
 * Primer:
 *   00000101 (5) -> 10100000 (160)
 */
static unsigned bitrev8(unsigned x) {
#pragma HLS INLINE
    unsigned r = 0;
    for (int i = 0; i < 8; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

/*
 * twiddle256
 * ----------
 * Vrne twiddle faktor W_256^k = exp(-j*2πk/256)
 * Uporabi COS_Q tabelo in simetrije
 * za izračun cos in -sin brez sin() funkcije.
 */
static cplx_t twiddle256(int k) {
#pragma HLS INLINE
    cplx_t w;

    if (k <= 64) {
        // 0..90°
        w.re = COS_Q[k];
        w.im = -COS_Q[64 - k];
    } else {
        // 90..180°
        int kp = 128 - k;
        w.re = -COS_Q[kp];
        w.im = -COS_Q[64 - kp];
    }
    return w;
}

/*
 * cadd
 * ----
 * Kompleksno seštevanje: (a + b)
 */
static cplx_t cadd(const cplx_t &a, const cplx_t &b) {
#pragma HLS INLINE
    return {a.re + b.re, a.im + b.im};
}

/*
 * csub
 * ----
 * Kompleksno odštevanje: (a - b)
 */
static cplx_t csub(const cplx_t &a, const cplx_t &b) {
#pragma HLS INLINE
    return {a.re - b.re, a.im - b.im};
}

/*
 * cmul
 * ----
 * Kompleksno množenje: a * b
 * Implementirano timing-friendly:
 *  - 4 DSP množenja
 *  - nato seštevanje/odštevanje
 *
 * (a.re + j a.im) * (b.re + j b.im)
 */
static void cmul(const cplx_t &a, const cplx_t &b, cplx_t &res) {
#pragma HLS INLINE off
#pragma HLS PIPELINE II=1
#pragma HLS LATENCY min=2 max=3

    fx_t ac = a.re * b.re;
    fx_t bd = a.im * b.im;
    fx_t ad = a.re * b.im;
    fx_t bc = a.im * b.re;

#pragma HLS BIND_OP variable=ac op=mul impl=dsp
#pragma HLS BIND_OP variable=bd op=mul impl=dsp
#pragma HLS BIND_OP variable=ad op=mul impl=dsp
#pragma HLS BIND_OP variable=bc op=mul impl=dsp

    res.re = ac - bd;
    res.im = ad + bc;
}

/*
 * axis_unpack
 * -----------
 * Razpakira AXI-Stream 32-bit besedo v kompleksno število.
 * Format:
 *   [15:0]  = real
 *   [31:16] = imag
 */
static cplx_t axis_unpack(const axis32_t &w) {
#pragma HLS INLINE
    ap_uint<32> d = w.data;

    ap_int<16> re_i = (ap_int<16>)d.range(15, 0);
    ap_int<16> im_i = (ap_int<16>)d.range(31, 16);

    cplx_t x;
    x.re.range() = re_i;
    x.im.range() = im_i;
    return x;
}

/*
 * axis_pack
 * ---------
 * Zapakira kompleksno število v AXI-Stream 32-bit format.
 * Signal last je 1 samo pri zadnjem vzorcu.
 */
static axis32_t axis_pack(const cplx_t &x, bool last) {
#pragma HLS INLINE
    ap_uint<32> d = 0;

    ap_int<16> re_i = (ap_int<16>)x.re.range();
    ap_int<16> im_i = (ap_int<16>)x.im.range();

    d.range(15, 0)  = (ap_uint<16>)re_i;
    d.range(31, 16) = (ap_uint<16>)im_i;

    axis32_t w;
    w.data = d;
    w.last = last ? 1 : 0;
    return w;
}

/*
 * read_stream_to_buf
 * ------------------
 * Prebere N kompleksnih vzorcev iz AXIS streama
 * in jih shrani v lokalni buffer.
 */
static void read_stream_to_buf(hls::stream<axis32_t> &s_in,
                               cplx_t buf[N]) {
#pragma HLS INLINE off
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        axis32_t w = s_in.read();
        buf[i] = axis_unpack(w);
    }
}

/*
 * write_buf_to_stream
 * -------------------
 * Zapiše N kompleksnih vzorcev v AXIS stream.
 * Zadnji vzorec ima last=1.
 */
static void write_buf_to_stream(const cplx_t buf[N],
                                hls::stream<axis32_t> &s_out) {
#pragma HLS INLINE off
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        s_out.write(axis_pack(buf[i], (i == N-1)));
    }
}

/*
 * fft_compute
 * -----------
 * Glavni FFT algoritem (radix-2, decimation-in-time).
 * 1) bit-reversal
 * 2) 8 FFT stopenj (za N=256)
 * 3) kopiranje rezultata
 */
static void fft_compute(const cplx_t inbuf[N],
                        cplx_t outbuf[N]) {
#pragma HLS INLINE off

    cplx_t buf[N];

    // Bit-reversal permutacija
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        buf[bitrev8((unsigned)i)] = inbuf[i];
    }

    // FFT stopnje
    for (int stage = 0; stage < 8; stage++) {
        int m       = 1 << (stage + 1);
        int half    = m >> 1;
        int tw_step = N / m;

        for (int grp = 0; grp < N; grp += m) {
            for (int j = 0; j < half; j++) {
#pragma HLS PIPELINE II=5

                int idxA = grp + j;
                int idxB = idxA + half;

                cplx_t u = buf[idxA];
                cplx_t w = twiddle256(j * tw_step);

                cplx_t t;
                cmul(buf[idxB], w, t);

                buf[idxA] = cadd(u, t);
                buf[idxB] = csub(u, t);
            }
        }
    }

    // Kopiranje rezultata
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        outbuf[i] = buf[i];
    }
}

/*
 * fft256
 * ------
 * Top-level HLS funkcija.
 * AXI-Stream vhod → FFT → AXI-Stream izhod
 */
void fft256(hls::stream<axis32_t> &s_in,
            hls::stream<axis32_t> &s_out) {
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE axis port=s_in
#pragma HLS INTERFACE axis port=s_out
#pragma HLS DATAFLOW

    cplx_t inbuf[N];
    cplx_t outbuf[N];

    read_stream_to_buf(s_in, inbuf);
    fft_compute(inbuf, outbuf);
    write_buf_to_stream(outbuf, s_out);
}


/*
===============================================================================
                               FFT256 BLOK DIAGRAM
===============================================================================

Top-level (DATAFLOW):
---------------------

   AXI4-Stream IN (axis32_t)
            |
            v
 +---------------------------+
 | read_stream_to_buf()      |
 |  - N=256 read             |
 |  - axis_unpack():         |
 |     [15:0]=Re, [31:16]=Im |
 +---------------------------+
            |
            v
 +---------------------------+
 | fft_compute()             |
 |  - bit-reversal reorder   |
 |  - 8 radix-2 FFT stages   |
 |  - output copy            |
 +---------------------------+
            |
            v
 +---------------------------+
 | write_buf_to_stream()     |
 |  - N=256 write            |
 |  - axis_pack():           |
 |     [15:0]=Re, [31:16]=Im |
 |  - last=1 on i==255       |
 +---------------------------+
            |
            v
   AXI4-Stream OUT (axis32_t)


===============================================================================
                        NOTRANJOST: fft_compute()
===============================================================================

1) Bit-reversal permutacija:
----------------------------

   inbuf[i]  --------------------->  buf[ bitrev8(i) ]

   (preuredi vzorce v bit-reversed vrstni red, da je radix-2 DIT FFT "prav"
    za in-place butterfly operacije)


2) 8 FFT stopenj (radix-2 DIT):
--------------------------------

Za stage = 0..7:

   m       = 2^(stage+1)          (velikost skupine/butterfly bloka)
   half    = m/2                  (razdalja med paroma A in B)
   tw_step = N/m                  (korak po twiddle faktorjih)

Znotraj vsake skupine (grp) in vsakega j (0..half-1):

         idxA = grp + j
         idxB = idxA + half

         u = buf[idxA]
         w = W_N^(j*tw_step)      (twiddle256())
         t = buf[idxB] * w        (cmul())

         buf[idxA] = u + t        (cadd)
         buf[idxB] = u - t        (csub)


ASCII “butterfly” za en par (idxA, idxB):
-----------------------------------------

                 u = buf[idxA]
                      |
                      |-----------------------> (+) -----> buf[idxA] = u + t
                      |                         ^
                      |                         |
                      |                       t = buf[idxB] * w
                      |                         |
                      v                         |
              buf[idxB] -----> (x w) -----------+
                      |
                      +-----------------------> (-) -----> buf[idxB] = u - t


3) Copy out:
------------

   outbuf[i] = buf[i]    za i=0..255


===============================================================================
                     TWIDDLE GENERACIJA (twiddle256)
===============================================================================

   twiddle256(k) vrne:  W_256^k = cos(2πk/256) - j*sin(2πk/256)

   COS_Q LUT ima cos za 0..90°; sin dobimo iz cos(90°-θ),
   ostale kvadrante pa s predznaki (simetrija).

===============================================================================
                     AXIS PAKIRANJE (32-bit)
===============================================================================

   axis32_t.data:
       bits [15:0]  -> Re (signed 16-bit)
       bits [31:16] -> Im (signed 16-bit)

   axis32_t.last:
       1 pri zadnjem vzorcu (i == 255), sicer 0

===============================================================================
*/