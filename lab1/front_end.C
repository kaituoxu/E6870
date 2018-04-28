
//  $Id: front_end.C,v 1.2 2016/01/23 03:15:23 stanchen Exp $

#include "front_end.H"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** Module for doing windowing. **/
void FrontEnd::do_window(const matrix<double>& in_feats,
                         matrix<double>& out_feats) const {
  //  Get parameters.
  //  Input samples per second.
  double sample_rate = get_float_param(params_, "window.sample_rate", 20000.0);
  //  Output frames per second.
  double frames_per_sec =
      get_float_param(params_, "window.frames_per_sec", 100.0);
  //  Width of each window, in seconds.
  double window_width = get_float_param(params_, "window.window_size", 0.025);
  //  Whether to do Hamming or rectangular windowing.
  bool do_Hamming = get_bool_param(params_, "window.hamming", true);

  //  Get number of input samples.
  int in_samp_cnt = in_feats.size1();
  if (in_feats.size2() != 1)
    throw runtime_error("Windowing expected vector input.");

  //  Input sampling period in seconds.
  double sample_period = 1.0 / sample_rate;
  //  Output frame period, in seconds.
  double frame_period = 1.0 / frames_per_sec;
  //  Number of samples per window.
  int samp_per_window = (int)(window_width / sample_period + 0.5);
  //  Number of samples to shift between each window.
  int samp_shift = (int)(frame_period / sample_period + 0.5);
  //  Number of output frames.
  int out_frame_cnt = (in_samp_cnt - samp_per_window) / samp_shift + 1;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(out_frame_cnt, samp_per_window);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      "in_feats", a matrix containing a single column holding the
  //      input samples for an utterance.  Each row holds a single sample.
  //
  //      in_feats(0 .. (in_samp_cnt - 1), 0)
  //
  //  Output:
  //      "out_feats", which should contain the result of windowing.
  //
  //      out_feats(0 .. (out_frame_cnt - 1), 0 .. (samp_per_window - 1))
  //
  //      Each row corresponds to a frame, and should hold the
  //      windowed samples for that frame.
  //      It has already been allocated to be of the correct size.
  //      If the boolean "do_Hamming" is true, then a Hamming
  //      window should be applied; otherwise, a rectangular
  //      window should be used.
  //
  //  See "in_samp_cnt", "samp_per_window", "samp_shift", and "out_frame_cnt"
  //  above for quantities you may (or may not) need for this computation.
  //
  //  When accessing matrices such as "in_feats" and "out_feats",
  //  use a syntax like "in_feats(frm_idx, dim_idx)" to access elements;
  //  using square brackets as in normal C arrays won't work.

  cout << format("samp_shift %d\n") % samp_shift;
  cout << format("out_frame_cnt %d\n") % out_frame_cnt;
  cout << format("samp_perWindw %d\n") % samp_per_window;
  cout << format("in_samp_cnt %d\n") %in_samp_cnt;
  if (do_Hamming) {
    // Hamming windows
    for (int r = 0; r < out_frame_cnt; ++r) {
      for (int c = 0; c < samp_per_window; ++c) {
        out_feats(r, c) =
            (0.54 - 0.46 * cos(2 * M_PI * c / (samp_per_window - 1))) *
            in_feats(r * samp_shift + c, 0);
      }
    }
  } else {
    // Rectangular window
    for (int r = 0; r < out_frame_cnt; ++r) {
      for (int c = 0; c < samp_per_window; ++c) {
        out_feats(r, c) = in_feats(r * samp_shift + c, 0);
      }
    }
  }
  //  END_LAB
}

/** Module for doing FFT. **/
void FrontEnd::do_fft(const matrix<double>& in_feats,
                      matrix<double>& out_feats) const {
  //  Make output dimension the smallest power of 2 at least as
  //  large as input dimension.
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = 2;
  while (out_dim_cnt < in_dim_cnt) out_dim_cnt *= 2;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  Input:
  //      "in_feats", a matrix with each row holding the windowed
  //      values for that frame.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      "out_feats", where an FFT should be applied to each
  //      row/frame of "in_feats".
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      For a given row/frame "frm_idx", the real and imaginary
  //      parts of the FFT value for frequency i/(out_dim_cnt*T)
  //      where T is the sample period are held in
  //      out_feats(frm_idx, 2*i) and out_feats(frm_idx, 2*i+1),
  //      respectively.

  vector<double> fft_buf;
  for (int frm_idx = 0; frm_idx < in_frame_cnt; ++frm_idx) {
    copy_matrix_row_to_vector(in_feats, frm_idx, fft_buf);
    //  Pad window with zeros, if needed.
    fft_buf.resize(out_dim_cnt, 0.0);
    real_fft(fft_buf);
    copy_vector_to_matrix_row(fft_buf, out_feats, frm_idx);
  }
}

/** Module for mel binning. **/
// change to google name style
void FrontEnd::do_melbin(const matrix<double>& in_feats,
                         matrix<double>& out_feats) const {
  //  Number of mel bins to make.
  int num_bins = get_int_param(params_, "melbin.bins", 26);
  //  Whether to take log of output or not.
  bool do_log = get_bool_param(params_, "melbin.log", true);
  //  Input samples per second.
  double sample_rate = get_float_param(params_, "window.sample_rate", 20000.0);
  double sample_period = 1.0 / sample_rate;

  //  Retrieve number of frames and dimension of input feature vectors.
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = num_bins;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      "in_feats", holding the output of a real FFT.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      "out_feats", which should contain the result of
  //      mel-binning.
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      If the boolean "doLog" is true,
  //      then each value should be replaced with its natural
  //      logarithm, or 0 if its logarithm is negative.
  //      "out_feats" has been allocated to be of the correct size.
  //
  //  See "in_frame_cnt", "in_dim_cnt", "out_dim_cnt", and "sample_period"
  //  above for quantities you will need for this computation.
  cout << format("sample_period %ls\n") % sample_period;
  cout << format("input dim %d\n") % in_dim_cnt;

  // TODO: Optimize the compuatation order to avoid duplicate computation
  int N = in_dim_cnt;
  int M = out_dim_cnt;
  double T = sample_period;
  // r means row
  for (int r = 0; r < in_frame_cnt; ++r) {
    // for each Mel bin, m is 1-based
    for (int m = 1; m <= M; ++m) {
      double sum = 0;
      // for each frequency (Hz)
      for (int i = 0; i < N / 2; ++i) {
        // BEGIN X(f)
        double f = i / (N * T);
        double real = in_feats(r, 2*i), img = in_feats(r, 2*i+1);  // S_2i, S_2i+1
        double X_f = sqrt(real*real + img*img);  // |X(f)|
        // END X(f)
        // BEGIN Hm[Mel_f]
        double Mel_f = 1127 * log(1 + f / 700);  // Mel(f)
        double Mel_f_max = 1127 * log(1 + 1 / (700 * 2 * T));
        double Mel_f_m = m * Mel_f_max / (M+1);  // Mel_f_min is 0
        double Mel_f_mp = (m-1) * Mel_f_max / (M+1); // p means previous
        double Mel_f_mn = (m+1) * Mel_f_max / (M+1); // n means next
        double H;
        if (Mel_f < Mel_f_mp || Mel_f > Mel_f_mn) {
          H = 0;
        } else if (Mel_f_mp <= Mel_f && Mel_f <= Mel_f_m) {
          H = (Mel_f - Mel_f_mp) / (Mel_f_m - Mel_f_mp);
        } else if (Mel_f_m <= Mel_f && Mel_f <= Mel_f_mn) {
          H = (Mel_f - Mel_f_mn) / (Mel_f_m - Mel_f_mn);
        } else {
          std::cout << "Invalid Mel(f) value!!" << std::endl;
        }
        // END Hm[Mel_f]
        sum += X_f * H;
      } // end for i
      // m is 1-based, but in out_feats, its position is 0-based
      if (do_log) {
        out_feats(r, m-1) = log(sum); // natrual logrithm
      } else {
        out_feats(r, m-1) = sum;
      }
    } // end for m
  } // end for r
  //  END_LAB
}

/** Module for doing discrete cosine transform. **/
// change to google name style
void FrontEnd::do_dct(const matrix<double>& in_feats,
                      matrix<double>& out_feats) const {
  //  Number of DCT coefficients to output.
  int num_coeffs = get_int_param(params_, "dct.coeffs", 12);
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = num_coeffs;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      The matrix "in_feats", holding the output of mel-binning.
  //
  //      in_feats(0 .. (in_frame_cnt - 1), 0 .. (in_dim_cnt - 1))
  //
  //  Output:
  //      The matrix "out_feats", which should contain the result of
  //      applying the DCT.
  //
  //      out_feats(0 .. (in_frame_cnt - 1), 0 .. (out_dim_cnt - 1))
  //
  //      "out_feats" has been allocated to be of the correct size.
  //
  //  See "in_frame_cnt", "in_dim_cnt", and "out_dim_cnt" above
  //  for quantities you will need for this computation.
  int N = in_dim_cnt;
  for (int r = 0; r < in_frame_cnt; ++r) {
    for (int j = 0; j < out_dim_cnt; ++j) {
      double sum = 0;
      for (int i = 0; i < N; ++i) {
        sum += in_feats(r, i) * cos(M_PI * (j + 1) * (i + 0.5) / N);
      }
      sum *= sqrt(2.0 / N);
      out_feats(r, j) = sum;
    }
  }
  //  END_LAB
}

/** Main signal processing routine.
 *   Calls each signal processing module in turn, unless
 *   parameter says not to.
 **/
void FrontEnd::get_feats(const matrix<double>& inAudio,
                         matrix<double>& out_feats) const {
  if (get_bool_param(params_, "frontend.null", false)) {
    out_feats = inAudio;
    return;
  }
  matrix<double> curFeats(inAudio);
  if (get_bool_param(params_, "frontend.window", true)) {
    do_window(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.fft", true)) {
    do_fft(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.melbin", true)) {
    do_melbin(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  if (get_bool_param(params_, "frontend.dct", true)) {
    do_dct(curFeats, out_feats);
    out_feats.swap(curFeats);
  }
  out_feats.swap(curFeats);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
