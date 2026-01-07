`timescale 1ns/1ps

module tb_fft256_sigwrap;

  localparam int N = 256;

  logic ap_clk = 0;
  always #5 ap_clk = ~ap_clk; // 100 MHz

  logic ap_rst_n;
  logic start;
  logic done;

  // -------------------------
  // INPUT TABLES (TB drives)
  // -------------------------
  logic [15:0] tb_in_re [0:N-1];
  logic [15:0] tb_in_im [0:N-1];

  fft256_sigwrap #(.N(N)) u_wrap (
    .ap_clk  (ap_clk),
    .ap_rst_n(ap_rst_n),
    .start   (start),
    .done    (done),

    // NEW: drive input tables into wrapper
    .in_re_q (tb_in_re),
    .in_im_q (tb_in_im)
  );

  // -------------------------
  // CSV paths 
  // -------------------------
  string path_time = "C:/Users/hp/Desktop/faks/DIVS/Projekt/plot/time_domain_sv.csv";
  string path_spec = "C:/Users/hp/Desktop/faks/DIVS/Projekt/plot/spectrum_sv.csv";

  integer f_time, f_spec;

  // -------------------------
  // Q8.8 helpers (match wrapper scaling)
  // -------------------------
function automatic logic [15:0] q8_8_from_real (real x);
  integer v;
  begin
    v = $rtoi(x * 256.0);
    if (v >  32767) v =  32767;
    if (v < -32768) v = -32768;
    q8_8_from_real = v[15:0];  
  end
endfunction



  function automatic real real_from_q8_8 (logic [15:0] q);
    integer sv;
    begin
      sv = $signed(q);
      real_from_q8_8 = sv / 256.0;
    end
  endfunction

  // -------------------------
  // Build 2-tone complex input ONCE
  // -------------------------
  initial begin : init_input
    int n;
    int k1, k2;
    real A1, A2;
    real dc_re, dc_im;
    real ph1, ph2;
    real re, im;

    // two frequencies (FFT bins)
k1 = 13;  A1 = 0.5;
k2 = 40;  A2 = 0.25;

for (n = 0; n < N; n++) begin
  ph1 = 2.0 * 3.141592653589793 * k1 * n / N;
  ph2 = 2.0 * 3.141592653589793 * k2 * n / N;

  re = A1*$cos(ph1) + A2*$cos(ph2);
  im = A1*$sin(ph1) + A2*$sin(ph2);        // IMAGINARNI DEL
 // im = 0.0;        //REALNI SIGNAL nima imaginarne komponente

  tb_in_re[n] = q8_8_from_real(re);
  tb_in_im[n] = q8_8_from_real(im);
end
  end

  task automatic start_run;
    begin
      @(posedge ap_clk);
      start <= 1'b1;
      repeat (20) @(posedge ap_clk);
      @(posedge ap_clk);
      start <= 1'b0;
    end
  endtask

  // -------------------------
  // MAIN
  // -------------------------
  initial begin
    int n;
    int k;
    real re, im, mag, mag_over_N, power;
    int bestk;
    real bestmag;

    // init
    ap_rst_n = 1'b0;
    start    = 1'b0;

    // short reset
    repeat (20) @(posedge ap_clk);
    @(posedge ap_clk);
    ap_rst_n <= 1'b1;
    $display("[%0t] Reset released", $time);

    repeat (20) @(posedge ap_clk);

    $display("[%0t] Starting run", $time);
    start_run();

    // wait done with timeout
    fork
      begin
        wait (done == 1'b1);
        $display("[%0t] DONE", $time);
      end
      begin
        #500us;
        $display("[%0t] TIMEOUT waiting for done", $time);
        $finish;
      end
    join_any
    disable fork;

    // -------------------------
    // Dump TIME domain CSV (INPUT)
    // -------------------------
    f_time = $fopen(path_time, "w");
    if (f_time == 0) begin
      $display("ERROR: cannot open %s", path_time);
      $finish;
    end

    $fwrite(f_time, "n,re,im\n");
    for (n = 0; n < N; n++) begin
      re = real_from_q8_8(tb_in_re[n]);
      im = real_from_q8_8(tb_in_im[n]);
      $fwrite(f_time, "%0d,%f,%f\n", n, re, im);
    end
    $fclose(f_time);

    // -------------------------
    // Dump SPECTRUM CSV (OUTPUT)
    // -------------------------
    f_spec = $fopen(path_spec, "w");
    if (f_spec == 0) begin
      $display("ERROR: cannot open %s", path_spec);
      $finish;
    end

    $fwrite(f_spec, "k,re,im,mag,mag_over_N,power\n");
    for (k = 0; k < N; k++) begin
      re = real_from_q8_8(u_wrap.out_re_q[k]);
      im = real_from_q8_8(u_wrap.out_im_q[k]);
      mag = $sqrt(re*re + im*im);
      mag_over_N = mag / N;
      power = re*re + im*im;
      $fwrite(f_spec, "%0d,%f,%f,%f,%f,%f\n", k, re, im, mag, mag_over_N, power);
    end
    $fclose(f_spec);

    // -------------------------
    // Find PEAK (rough check)
    // -------------------------
    bestk = 0;
    bestmag = -1.0;
    for (k = 0; k < N; k++) begin
      re = real_from_q8_8(u_wrap.out_re_q[k]);
      im = real_from_q8_8(u_wrap.out_im_q[k]);
      mag = $sqrt(re*re + im*im);
      if (mag > bestmag) begin
        bestmag = mag;
        bestk = k;
      end
    end

    $display("PEAK bin=%0d  |X|/N=%f", bestk, bestmag / N);
    $display("Wrote CSV:\n  %s\n  %s", path_time, path_spec);

    repeat (50) @(posedge ap_clk);
    $finish;
  end

endmodule
