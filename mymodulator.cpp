#ifndef MYMODULATOR
#include "mymodulator.h"
#endif

void MyQAM::set_M(int Mary)
{
  k = levels2bits(Mary);
  M = Mary;
  it_assert((pow2i(k) == M) && (is_even(k)),
            "QAM::set_M(): M = " << M << " is not an even power of 2");
  L = round_i(std::sqrt(static_cast<double>(M)));

  double average_energy = (M - 1) * 2.0 / 3.0;
  scaling_factor = std::sqrt(average_energy);

  symbols.set_size(M);
  bitmap.set_size(M, k);
  bits2symbols.set_size(M);

  bmat gray_code = graycode(levels2bits(L));


  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      symbols(i*L + j) = std::complex<double>(((L - 1) - j * 2) / scaling_factor,
                                              ((L - 1) - i * 2) / scaling_factor);
      bitmap.set_row(i*L + j, concat(gray_code.get_row(i),
                                     gray_code.get_row(j)));
      bits2symbols(bin2dec(bitmap.get_row(i*L + j))) = i * L + j;
    }
  }
  calculate_softbit_matrices();

  setup_done = true;
}



void MyQAM::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "QAM::demodulate_bits(): Modulator not ready.");
  out.set_size(k*signal.size(), false);

  int temp_real, temp_imag;

  for (int i = 0; i < signal.size(); i++) {
    temp_real = round_i((L - 1) - (std::real(signal(i) * scaling_factor)
                                   + (L - 1)) / 2.0);
    temp_imag = round_i((L - 1) - (std::imag(signal(i) * scaling_factor)
                                   + (L - 1)) / 2.0);
    if (temp_real < 0)
      temp_real = 0;
    else if (temp_real > (L - 1))
      temp_real = (L - 1);
    if (temp_imag < 0)
      temp_imag = 0;
    else if (temp_imag > (L - 1))
      temp_imag = (L - 1);
    out.replace_mid(k*i, bitmap.get_row(temp_imag * L + temp_real));
  }
}

bvec MyQAM::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}



void Paper19QAM::set_M()
{
  //k = levels2bits(Mary);
  k = 4;
  M = 16;//Mary;
  /*
  it_assert((pow2i(k) == M) && (is_even(k)),
            "Paper19QAM::set_M(): M = " << M << " is not an even power of 2");
  */
  L = round_i(std::sqrt(static_cast<double>(M)));
  symbols.set_size(M);
  bitmap.set_size(M, k);
  bits2symbols.set_size(M);

  // symbols：symbols(符号)=信号
  symbols(7) = std::complex<double>(0,0);
  symbols(8) = std::complex<double>(1,0);
  symbols(3) = std::complex<double>(-0.5, sqrt(3)/2) * symbols(8);
  symbols(4) = std::complex<double>(0.5, sqrt(3)/2) * symbols(8);
  symbols(6) = -symbols(8);
  symbols(11) = -symbols(4);
  symbols(12) = -symbols(3);

  symbols(0) = 2 * symbols(3);
  symbols(9) = 2 * symbols(8);
  symbols(14) = 2 * symbols(11);

  symbols(15) = symbols(11)+symbols(12);
  symbols(2) = symbols(3) + symbols(6);
  symbols(5) = symbols(8) + symbols(4);
  symbols(10) = symbols(11) + symbols(6);
  symbols(13) = symbols(8) + symbols(12);
  symbols(1) = symbols(3) + symbols(4);

  // bitmap：bitmap.set_row(符号)=bit列
  bitmap.set_row(0,bvec("1 0 1 0"));
  bitmap.set_row(1,bvec("0 1 1 0"));
  bitmap.set_row(2,bvec("1 1 1 1"));
  bitmap.set_row(3,bvec("0 0 1 1"));
  bitmap.set_row(4,bvec("0 0 1 0"));
  bitmap.set_row(5,bvec("0 1 1 1"));
  bitmap.set_row(6,bvec("0 1 0 1"));
  bitmap.set_row(7,bvec("0 0 0 0"));
  bitmap.set_row(8,bvec("0 0 0 1"));
  bitmap.set_row(9,bvec("1 0 1 1"));
  bitmap.set_row(10,bvec("1 1 0 0"));
  bitmap.set_row(11,bvec("0 1 0 0"));
  bitmap.set_row(12,bvec("1 0 0 0"));
  bitmap.set_row(13,bvec("1 0 0 1"));
  bitmap.set_row(14,bvec("1 1 0 1"));
  bitmap.set_row(15,bvec("1 1 1 0"));

  double sum_energy = 0;
  for (int i = 0; i < M; i++) {
    // bits2symbols：bits2symbols(bit列の10進数)=符号
    // 基本は上のsymbolsとbitmapを対応づけているだけだから
    // 何もしなくてOK
    bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    sum_energy += std::pow(std::abs(symbols(i)),2);
  }

  double average_energy = sum_energy / M;
  scaling_factor = std::sqrt(average_energy);
  symbols /= scaling_factor;

  calculate_softbit_matrices();

  setup_done = true;
}

// MLD
void Paper19QAM::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "Paper19QAM::demodulate_bits(): Modulator not ready.");
  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
    int c_ind = 0;
    double dif = std::abs(signal(i) - symbols(c_ind));
    double tmp_dif;

    for (int c = 1; c < symbols.size(); c++) {
      tmp_dif = std::abs(signal(i) - symbols(c));
      if (tmp_dif < dif) {
        c_ind = c;
        dif = tmp_dif;
      }
    }

    out.replace_mid(k*i, bitmap.get_row(c_ind));
  }
}

bvec Paper19QAM::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}


void TriSq19QAM::set_M()
{
  //k = levels2bits(Mary);
  k = 4;
  M = 16;//Mary;
  /*
  it_assert((pow2i(k) == M) && (is_even(k)),
            "Paper19QAM::set_M(): M = " << M << " is not an even power of 2");
  */
  L = round_i(std::sqrt(static_cast<double>(M)));
  symbols.set_size(M);
  bitmap.set_size(M, k);
  bits2symbols.set_size(M);

  double A = 0.4;
  std::complex<double> L1 = std::complex<double>(-A / 2, sqrt(3) * A / 2);
  std::complex<double> R1 = std::complex<double>(A / 2, sqrt(3) * A / 2);
  std::complex<double> L2 = std::complex<double>(-sqrt(3) * A / 2, A / 2);
  std::complex<double> R2 = std::complex<double>(sqrt(3) * A / 2, A / 2);

  // symbols：symbols(符号)=信号
  std::complex<double> tmp = std::complex<double>(0, A / 2);

  symbols(0) = tmp + L1;
  symbols(1) = tmp + R1;
  symbols(2) = tmp - R2 + L1;
  symbols(3) = tmp - L2 + R1;
//  symbols(6) = tmp;
  symbols(6) = symbols(2) - R1;
  symbols(7) = tmp - R2;
  symbols(8) = tmp - L2;
  symbols(9) = tmp -L2 + R1 - L1;
  symbols(12) = tmp - R2 - R1;
  symbols(13) = tmp - L2 - L1;
  symbols(14) = -tmp - R1;
  symbols(15) = -tmp - L1;
  symbols(4) = (symbols(0) + symbols(2) + symbols(7) + tmp) / 4;
  symbols(5) = (symbols(1) + symbols(3) + symbols(8) + tmp) / 4;
  symbols(10) = (symbols(7) + symbols(12) + symbols(14) - tmp) / 4;
  symbols(11) = (symbols(8) + symbols(13) + symbols(15) - tmp) / 4;

  // bitmap：bitmap.set_row(符号)=bit列
  bitmap.set_row(0,bvec("1 0 1 0"));
  bitmap.set_row(1,bvec("1 1 0 0"));
  bitmap.set_row(2,bvec("0 0 1 0"));
  bitmap.set_row(3,bvec("0 1 0 1"));
  bitmap.set_row(4,bvec("0 0 1 1"));
  bitmap.set_row(5,bvec("0 1 0 0"));
  bitmap.set_row(6,bvec("1 1 1 0"));
  bitmap.set_row(7,bvec("1 1 1 1"));
  bitmap.set_row(8,bvec("0 0 0 0"));
  bitmap.set_row(9,bvec("0 0 0 1"));
  bitmap.set_row(10,bvec("0 1 1 1"));
  bitmap.set_row(11,bvec("1 0 0 0"));
  bitmap.set_row(12,bvec("0 1 1 0"));
  bitmap.set_row(13,bvec("1 0 0 1"));
  bitmap.set_row(14,bvec("1 0 1 1"));
  bitmap.set_row(15,bvec("1 1 0 1"));

  double sum_energy = 0;
  for (int i = 0; i < M; i++) {
    // bits2symbols：bits2symbols(bit列の10進数)=符号
    // 基本は上のsymbolsとbitmapを対応づけているだけだから
    // 何もしなくてOK
    bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    sum_energy += std::pow(std::abs(symbols(i)),2);
  }

  double average_energy = sum_energy / M;
  scaling_factor = std::sqrt(average_energy);
  symbols /= scaling_factor;

  calculate_softbit_matrices();

  setup_done = true;
}

// MLD
void TriSq19QAM::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "Paper19QAM::demodulate_bits(): Modulator not ready.");
  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
    int c_ind = 0;
    double dif = std::abs(signal(i) - symbols(c_ind));
    double tmp_dif;

    for (int c = 1; c < symbols.size(); c++) {
      tmp_dif = std::abs(signal(i) - symbols(c));
      if (tmp_dif < dif) {
        c_ind = c;
        dif = tmp_dif;
      }
    }

    out.replace_mid(k*i, bitmap.get_row(c_ind));
  }
}

bvec TriSq19QAM::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}
