`timescale 1ns/1ps

module fft256_sigwrap #(
  parameter int N = 256
)(
  input  logic ap_clk,
  input  logic ap_rst_n,
  input  logic start,       // level start (TB ga drži 1 nekaj ciklov)
  output logic done,

  // input tables come from outside (TB or other module) ----
  input  logic [15:0] in_re_q [0:N-1],
  input  logic [15:0] in_im_q [0:N-1]
);

  // ---------------- ap_ctrl_hs ----------------
  logic ap_start;
  logic ap_done, ap_idle, ap_ready;

  // ---------------- AXIS in ----------------
  logic [31:0] s_in_tdata;
  logic [3:0]  s_in_tkeep;
  logic [3:0]  s_in_tstrb;
  logic        s_in_tvalid;
  logic        s_in_tready;
  logic        s_in_tlast;

  // ---------------- AXIS out ----------------
  logic [31:0] s_out_tdata;
  logic [3:0]  s_out_tkeep;
  logic [3:0]  s_out_tstrb;
  logic        s_out_tvalid;
  logic        s_out_tready;
  logic        s_out_tlast;

  // ---------------- DUT ----------------
  fft256_1 dut (
    .ap_clk(ap_clk),
    .ap_rst_n(ap_rst_n),
    .ap_done(ap_done),
    .ap_idle(ap_idle),
    .ap_ready(ap_ready),
    .ap_start(ap_start),

    .s_in_TDATA (s_in_tdata),
    .s_in_TKEEP (s_in_tkeep),
    .s_in_TLAST (s_in_tlast),
    .s_in_TREADY(s_in_tready),
    .s_in_TSTRB (s_in_tstrb),
    .s_in_TVALID(s_in_tvalid),

    .s_out_TDATA (s_out_tdata),
    .s_out_TKEEP (s_out_tkeep),
    .s_out_TLAST (s_out_tlast),
    .s_out_TREADY(s_out_tready),
    .s_out_TSTRB (s_out_tstrb),
    .s_out_TVALID(s_out_tvalid)
  );

 
  logic [15:0] out_re_q[0:N-1];
  logic [15:0] out_im_q[0:N-1];

  int in_idx, out_idx;

  // ---------------- SIM-ONLY "Analog" probes ----------------
`ifndef SYNTHESIS
  real in_re_analog, in_im_analog;
  real out_re_analog, out_im_analog;

  function automatic real real_from_q8_8 (logic [15:0] q);
    integer sv;
    begin
      sv = $signed(q);
      real_from_q8_8 = sv / 256.0;
    end
  endfunction
`endif

  // ---------------- Start pulse detect ----------------
  logic start_d;
  logic start_pulse;

  always_ff @(posedge ap_clk) begin
    if (!ap_rst_n) start_d <= 1'b0;
    else          start_d <= start;
  end
  assign start_pulse = start & ~start_d;

  // ---------------- Constants ----------------
  always_comb begin
    s_in_tkeep = 4'hF;
    s_in_tstrb = 4'hF;
  end

  // Receive ready only during RECV
  typedef enum logic [1:0] {S_IDLE, S_SEND, S_RECV, S_DONE} state_t;
  state_t st;

  always_comb begin
    s_out_tready = (st == S_RECV);
  end

  // ---------------- Main FSM ----------------
  always_ff @(posedge ap_clk) begin
    if (!ap_rst_n) begin
      st <= S_IDLE;
      done <= 1'b0;

      // KEY: ap_start low in reset
      ap_start <= 1'b0;

      s_in_tvalid <= 1'b0;
      s_in_tlast  <= 1'b0;
      s_in_tdata  <= '0;

      in_idx <= 0;
      out_idx <= 0;

`ifndef SYNTHESIS
      in_re_analog  <= 0.0;
      in_im_analog  <= 0.0;
      out_re_analog <= 0.0;
      out_im_analog <= 0.0;
`endif

    end else begin
      done <= 1'b0;

      // KEY: Keep core enabled all the time (streaming-style use)
      ap_start <= 1'b1;

      case (st)
        S_IDLE: begin
          s_in_tvalid <= 1'b0;
          s_in_tlast  <= 1'b0;
          in_idx <= 0;
          out_idx <= 0;

          if (start_pulse) begin
            st <= S_SEND;
          end
        end

        S_SEND: begin
          // present a word when not already holding valid
          if (!s_in_tvalid && (in_idx < N)) begin
            s_in_tdata  <= {in_im_q[in_idx], in_re_q[in_idx]};
            s_in_tlast  <= (in_idx == (N-1));
            s_in_tvalid <= 1'b1;

`ifndef SYNTHESIS
            // analog probes update per sample (SIM only)
            in_re_analog <= real_from_q8_8(in_re_q[in_idx]);
            in_im_analog <= real_from_q8_8(in_im_q[in_idx]);
`endif
          end

          // handshake
          if (s_in_tvalid && s_in_tready) begin
            s_in_tvalid <= 1'b0;
            s_in_tlast  <= 1'b0;
            in_idx <= in_idx + 1;

            if (in_idx == (N-1)) begin
              st <= S_RECV;
            end
          end
        end

        S_RECV: begin
          if (s_out_tvalid && s_out_tready) begin
            out_re_q[out_idx] <= s_out_tdata[15:0];
            out_im_q[out_idx] <= s_out_tdata[31:16];

`ifndef SYNTHESIS
            out_re_analog <= real_from_q8_8(s_out_tdata[15:0]);
            out_im_analog <= real_from_q8_8(s_out_tdata[31:16]);
`endif

            out_idx <= out_idx + 1;

            if (out_idx == (N-1)) begin
              st <= S_DONE;
            end
          end
        end

        S_DONE: begin
          done <= 1'b1;
          st <= S_IDLE;
        end
      endcase
    end
  end

endmodule
