#include <cmath>
#include <cstdio>
#include <fstream>
#include "fft.h"

static cplx_t make_cplx(double re, double im)
{
    cplx_t x;
    x.re = (fx_t)re;
    x.im = (fx_t)im;
    return x;
}

// ---------------- AXIS pack / unpack ----------------
// AXIS payload: [15:0]=re(raw bits), [31:16]=im(raw bits)
static axis32_t pack_axis(const cplx_t &x, bool last)
{
    ap_uint<32> d = 0;
    ap_int<16> re_i = (ap_int<16>)x.re.range();
    ap_int<16> im_i = (ap_int<16>)x.im.range();

    d.range(15, 0)  = (ap_uint<16>)re_i;
    d.range(31,16)  = (ap_uint<16>)im_i;

    axis32_t w;
    w.data = d;
    w.last = last ? 1 : 0;
    return w;
}

static cplx_t unpack_axis(const axis32_t &w)
{
    ap_uint<32> d = w.data;
    ap_int<16> re_i = (ap_int<16>)d.range(15,0);
    ap_int<16> im_i = (ap_int<16>)d.range(31,16);

    cplx_t x;
    x.re.range() = re_i;
    x.im.range() = im_i;
    return x;
}
// ---------------------------------------------------

int main()
{
    hls::stream<axis32_t> s_in, s_out;

    cplx_t in[N], out[N];

    // ===================== SIGNAL (3 KOMPONENTE) =====================
    const int    k1 = 13;   const double A1 = 0.5;
    const int    k2 = 40;   const double A2 = 1.0;
    const int    k3 = 90;   const double A3 = 0.3;

    const double dc_re = 0.0;
    const double dc_im = 0.0;

    for (int n = 0; n < N; n++) {
        double re = dc_re;
       double im = dc_im;
       

        double ph1 = 2.0 * M_PI * (double)k1 * (double)n / (double)N;
        re += A1 * std::cos(ph1);
        //im += A1 * std::sin(ph1);    da je realen signal -------------------------------------
        im=0;


        double ph2 = 2.0 * M_PI * (double)k2 * (double)n / (double)N;
        re += A2 * std::cos(ph2);
       // im += A2 * std::sin(ph2);     da je realen signal ------------------
        im=0;
        double ph3 = 2.0 * M_PI * (double)k3 * (double)n / (double)n / (double)N; 

     
        ph3 = 2.0 * M_PI * (double)k3 * (double)n / (double)N;
        re += A3 * std::cos(ph3);
       // im += A3 * std::sin(ph3);          da je realen signal --------------
       im=0;           

        in[n] = make_cplx(re, im);

        // send to stream
        s_in.write(pack_axis(in[n], (n == N-1)));
    }

    // ===================== RUN FFT ====================
    fft256(s_in, s_out);

    // ===================== READ OUT ===================
    for (int i = 0; i < N; i++) {
        axis32_t w = s_out.read();
        out[i] = unpack_axis(w);
    }

    // ===================== CSV PATHS ==================
    const char* path_time = R"(C:\Users\hp\Desktop\faks\DIVS\Projekt\plot\time_domain.csv)";
    const char* path_spec = R"(C:\Users\hp\Desktop\faks\DIVS\Projekt\plot\spectrum.csv)";

    // ===================== CSV: TIME ==================
    {
        std::ofstream f(path_time);
        if (!f.is_open()) {
            printf("ERROR: cannot open %s\n", path_time);
            return 1;
        }
        f << "n,re,im\n";
        for (int n = 0; n < N; n++) {
            f << n << ","
              << (double)in[n].re << ","
              << (double)in[n].im << "\n";
        }
    }

    // ===================== CSV: FFT ===================
    {
        std::ofstream f(path_spec);
        if (!f.is_open()) {
            printf("ERROR: cannot open %s\n", path_spec);
            return 1;
        }

        f << "k,re,im,mag,mag_over_N,power\n";
        for (int k = 0; k < N; k++) {
            double re = (double)out[k].re;
            double im = (double)out[k].im;

            double mag = std::sqrt(re*re + im*im);
            double mag_over_N = mag / (double)N;
            double power = re*re + im*im;

            f << k << ","
              << re << ","
              << im << ","
              << mag << ","
              << mag_over_N << ","
              << power << "\n";
        }
    }

    // ===================== QUICK CHECK (amplitude) =================
    int bestk = 0;
    double bestmag = -1.0;

    for (int k = 0; k < N; k++) {
        double re = (double)out[k].re;
        double im = (double)out[k].im;
        double mag = std::sqrt(re*re + im*im);
        if (mag > bestmag) { bestmag = mag; bestk = k; }
    }

    printf("FFT peak bin = %d, |X|/N = %.6f\n", bestk, bestmag / (double)N);
    printf("Wrote:\n  %s\n  %s\n", path_time, path_spec);

    return 0;
}
