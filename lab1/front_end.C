
//  $Id: front_end.C,v 1.2 2016/01/23 03:15:23 stanchen Exp $

#include "front_end.H"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** Module for doing windowing. **/
void FrontEnd::do_window(const matrix<double>& inFeats,
                         matrix<double>& outFeats) const {
  //  Get parameters.
  //  Input samples per second.
  double sampleRate = get_float_param(m_params, "window.sample_rate", 20000.0);
  //  Output frames per second.
  double framesPerSec =
      get_float_param(m_params, "window.frames_per_sec", 100.0);
  //  Width of each window, in seconds.
  double windowWidth = get_float_param(m_params, "window.window_size", 0.025);
  //  Whether to do Hamming or rectangular windowing.
  bool doHamming = get_bool_param(m_params, "window.hamming", true);

  //  Get number of input samples.
  int inSampCnt = inFeats.size1();
  if (inFeats.size2() != 1)
    throw runtime_error("Windowing expected vector input.");

  //  Input sampling period in seconds.
  double samplePeriod = 1.0 / sampleRate;
  //  Output frame period, in seconds.
  double framePeriod = 1.0 / framesPerSec;
  //  Number of samples per window.
  int sampPerWindow = (int)(windowWidth / samplePeriod + 0.5);
  //  Number of samples to shift between each window.
  int sampShift = (int)(framePeriod / samplePeriod + 0.5);
  //  Number of output frames.
  int outFrameCnt = (inSampCnt - sampPerWindow) / sampShift + 1;

  //  Allocate output matrix and fill with zeros.
  outFeats.resize(outFrameCnt, sampPerWindow);
  outFeats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      "inFeats", a matrix containing a single column holding the
  //      input samples for an utterance.  Each row holds a single sample.
  //
  //      inFeats(0 .. (inSampCnt - 1), 0)
  //
  //  Output:
  //      "outFeats", which should contain the result of windowing.
  //
  //      outFeats(0 .. (outFrameCnt - 1), 0 .. (sampPerWindow - 1))
  //
  //      Each row corresponds to a frame, and should hold the
  //      windowed samples for that frame.
  //      It has already been allocated to be of the correct size.
  //      If the boolean "doHamming" is true, then a Hamming
  //      window should be applied; otherwise, a rectangular
  //      window should be used.
  //
  //  See "inSampCnt", "sampPerWindow", "sampShift", and "outFrameCnt"
  //  above for quantities you may (or may not) need for this computation.
  //
  //  When accessing matrices such as "inFeats" and "outFeats",
  //  use a syntax like "inFeats(frmIdx, dimIdx)" to access elements;
  //  using square brackets as in normal C arrays won't work.

  cout << format("sampShift %d\n") % sampShift;
  cout << format("outFrameCnt %d\n") % outFrameCnt;
  cout << format("sampPerWindw %d\n") % sampPerWindow;
  cout << format("inSampCnt %d\n") %inSampCnt;
  if (doHamming) {
    // Hamming windows
    for (int r = 0; r < outFrameCnt; ++r) {
      for (int c = 0; c < sampPerWindow; ++c) {
        outFeats(r, c) =
            (0.54 - 0.46 * cos(2 * M_PI * c / (sampPerWindow - 1))) *
            inFeats(r * sampShift + c, 0);
      }
    }
  } else {
    // Rectangular window
    for (int r = 0; r < outFrameCnt; ++r) {
      for (int c = 0; c < sampPerWindow; ++c) {
        outFeats(r, c) = inFeats(r * sampShift + c, 0);
      }
    }
  }
  //  END_LAB
}

/** Module for doing FFT. **/
void FrontEnd::do_fft(const matrix<double>& inFeats,
                      matrix<double>& outFeats) const {
  //  Make output dimension the smallest power of 2 at least as
  //  large as input dimension.
  int inFrameCnt = inFeats.size1();
  int inDimCnt = inFeats.size2();
  int outDimCnt = 2;
  while (outDimCnt < inDimCnt) outDimCnt *= 2;

  //  Allocate output matrix and fill with zeros.
  outFeats.resize(inFrameCnt, outDimCnt);
  outFeats.clear();

  //  Input:
  //      "inFeats", a matrix with each row holding the windowed
  //      values for that frame.
  //
  //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
  //
  //  Output:
  //      "outFeats", where an FFT should be applied to each
  //      row/frame of "inFeats".
  //
  //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
  //
  //      For a given row/frame "frmIdx", the real and imaginary
  //      parts of the FFT value for frequency i/(outDimCnt*T)
  //      where T is the sample period are held in
  //      outFeats(frmIdx, 2*i) and outFeats(frmIdx, 2*i+1),
  //      respectively.

  vector<double> fftBuf;
  for (int frmIdx = 0; frmIdx < inFrameCnt; ++frmIdx) {
    copy_matrix_row_to_vector(inFeats, frmIdx, fftBuf);
    //  Pad window with zeros, if needed.
    fftBuf.resize(outDimCnt, 0.0);
    real_fft(fftBuf);
    copy_vector_to_matrix_row(fftBuf, outFeats, frmIdx);
  }
}

/** Module for mel binning. **/
// change to google name style
void FrontEnd::do_melbin(const matrix<double>& in_feats,
                         matrix<double>& out_feats) const {
  //  Number of mel bins to make.
  int num_bins = get_int_param(m_params, "melbin.bins", 26);
  //  Whether to take log of output or not.
  bool do_log = get_bool_param(m_params, "melbin.log", true);
  //  Input samples per second.
  double sample_rate = get_float_param(m_params, "window.sample_rate", 20000.0);
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
  //      "inFeats", holding the output of a real FFT.
  //
  //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
  //
  //  Output:
  //      "outFeats", which should contain the result of
  //      mel-binning.
  //
  //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
  //
  //      If the boolean "doLog" is true,
  //      then each value should be replaced with its natural
  //      logarithm, or 0 if its logarithm is negative.
  //      "outFeats" has been allocated to be of the correct size.
  //
  //  See "inFrameCnt", "inDimCnt", "outDimCnt", and "samplePeriod"
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
  int num_coeffs = get_int_param(m_params, "dct.coeffs", 12);
  int in_frame_cnt = in_feats.size1();
  int in_dim_cnt = in_feats.size2();
  int out_dim_cnt = num_coeffs;

  //  Allocate output matrix and fill with zeros.
  out_feats.resize(in_frame_cnt, out_dim_cnt);
  out_feats.clear();

  //  BEGIN_LAB
  //
  //  Input:
  //      The matrix "inFeats", holding the output of mel-binning.
  //
  //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
  //
  //  Output:
  //      The matrix "outFeats", which should contain the result of
  //      applying the DCT.
  //
  //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
  //
  //      "outFeats" has been allocated to be of the correct size.
  //
  //  See "inFrameCnt", "inDimCnt", and "outDimCnt" above
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
                         matrix<double>& outFeats) const {
  if (get_bool_param(m_params, "frontend.null", false)) {
    outFeats = inAudio;
    return;
  }
  matrix<double> curFeats(inAudio);
  if (get_bool_param(m_params, "frontend.window", true)) {
    do_window(curFeats, outFeats);
    outFeats.swap(curFeats);
  }
  if (get_bool_param(m_params, "frontend.fft", true)) {
    do_fft(curFeats, outFeats);
    outFeats.swap(curFeats);
  }
  if (get_bool_param(m_params, "frontend.melbin", true)) {
    do_melbin(curFeats, outFeats);
    outFeats.swap(curFeats);
  }
  if (get_bool_param(m_params, "frontend.dct", true)) {
    do_dct(curFeats, outFeats);
    outFeats.swap(curFeats);
  }
  outFeats.swap(curFeats);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
