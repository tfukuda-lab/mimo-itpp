#include <itpp/itcomm.h>
#include <itpp/comm/modulator.h>
#ifndef MYMODULATOR
#include "mymodulator.h"
#endif

using namespace itpp;
//These lines are needed for use of cout and endl
using std::cout;
using std::endl;

int main()
{
  //Declarations of scalars and vectors:
  int i, Number_of_bits, Number_of_frames;
  int Nt = 20;
  int Nr = 20;
  double EbN0_lowest = 0.0;
  double EbN0_highest = 50;
  int EbN0_num = 11;
  cmat H,Y,X,N;
  double Ec, Eb;
  vec EbN0dB, EbN0, N0, noise_variance; //vec is a vector containing double
  bvec transmitted_bits, received_bits;                 //bvec is a vector containing bits
  cvec transmitted_symbols, received_symbols;           //cvec is a vector containing double_complex
  //Declarations of classes:
  Real_Timer tt;                 //The timer used to measure the execution time
  BERC berc;                     //Used to count the bit errors
  
  // to obtain BER
  int numSimulate = 3;
  MyQAM *qams[numSimulate] = {new MyQAM(16), new Paper19QAM(), new TriSq19QAM()};
  vec bit_error_rates[numSimulate];                     //Used to count the bit errors
  std::string labels[numSimulate] = {"16QAM", "Paper19QAM", "TriSquare"};

  //Reset and start the timer:
  tt.tic();

  //Init:
  Ec = 1.0;                      //The transmitted energy per QPSK symbol is 1.
  Eb = Ec / 4.0;                 //The transmitted energy per bit is 0.5.
  EbN0dB = linspace(EbN0_lowest, EbN0_highest, EbN0_num); //Simulate for 10 Eb/N0 values from 0 to 9 dB.
  EbN0 = inv_dB(EbN0dB);         //Calculate Eb/N0 in a linear scale instead of dB.
  N0 = Eb * pow(EbN0, -1.0);     //N0 is the variance of the (complex valued) noise.
  Number_of_bits = 10000;       //One hundred thousand bits is transmitted for each Eb/N0 value
  Number_of_frames = 10000;       //One hundred thousand bits is transmitted for each Eb/N0 value

  for (int i = 0; i < numSimulate; i++) {
    bit_error_rates[i].set_size(EbN0dB.length(), false);
  }
  //Randomize the random number generators in it++:
  RNG_randomize();

  //Iterate over all EbN0dB values:
  for (i = 0; i < EbN0dB.length(); i++) {
    //Show how the simulation progresses:
    cout << "Now simulating Eb/N0 value number " << i + 1 << " of " << EbN0dB.length() <<"(var=" << N0(i) << ")" << endl;

    //Generate a vector of random bits to transmit:
    transmitted_bits = randb(Number_of_bits);

    for (int j = 0; j < numSimulate; j++) {
        //Modulate the bits to 16QAM symbols:
        transmitted_symbols = qams[j]->modulate_bits(transmitted_bits);
        int n_timeslot = transmitted_symbols.size() / Nt;
        X = reshape(transmitted_symbols, Nt, n_timeslot);
        berc.clear();                               //Clear the bit error rate counter
        for (int f = 0; f < Number_of_frames; f++) {
            H = randn_c(Nr, Nt) / sqrt(2);
            N = randn_c(Nr, n_timeslot) * std::complex<double>(sqrt(N0(i) / 2));// * std::complex<double>(0);
            Y = H * X + N;

            cmat Hh = H.hermitian_transpose();
            cmat Wzf = inv(Hh * H)*Hh;
            cmat Z = Wzf * Y;

            received_symbols = rvectorize(Z.transpose());
            //Demodulate the received QPSK symbols into received bits:
            received_bits = qams[j]->demodulate_bits(received_symbols);
            //Calculate the bit error rate:
            berc.count(transmitted_bits, received_bits); //Count the bit errors
        }
        bit_error_rates[j](i) = berc.get_errorrate();   //Save the estimated BER in the result vector
    }
  }

  tt.toc();

  // show the result
  for (int i = 0; i < numSimulate; i++) {
    //Print the results:
    cout << endl;
    cout << "EbN0dB = " << EbN0dB << " [dB]" << endl;
    cout << "BER = " << bit_error_rates[i] << endl;
    cout << "Saving results to ./" << labels[i] << "_result_file.it" << endl;
    cout << endl;
  }

  // save the graph

  FILE *gp;
  gp = popen("gnuplot", "w");
  //fprintf(gp, "unset key\n");
  fprintf(gp, "set terminal png\n");
  fprintf(gp, "set output 'result.png'\n");

  fprintf(gp, "set logscale y\n");
  //fprintf(gp, "set xrange[%lf:%lf]\n",EbN0_lowest,EbN0_highest);
  //fprintf(gp, "set yrange[0.00000001:1]\n");
  fprintf(gp, "set grid\n");
  fprintf(gp, "set xlabel \"EbN0\"\n");
  fprintf(gp, "set ylabel \"BER\"\n");

  fprintf(gp, "plot '-' with lines linetype 1 title '16QAM', '-' with lines linetype 2 title '6', '-' with lines linetype 3 title '3-4'\n");

  for (int i = 0; i < numSimulate; i++) {
    for (int x = 0; x < EbN0dB.size(); x++) {
      fprintf(gp, "%lf\t%lf\n", EbN0dB(x), bit_error_rates[i](x));
    }
    fprintf(gp, "e\n");
  }

  fflush(gp);

  //Exit program:
  return 0;
}