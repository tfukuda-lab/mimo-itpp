#define MYMODULATOR

#include <itpp/itcomm.h>
#include <itpp/comm/modulator.h>

using namespace itpp;
//These lines are needed for use of cout and endl
using std::cout;
using std::endl;



class ITPP_EXPORT MyQAM : public Modulator<std::complex<double> >
{
public:
  //! Default Constructor
  MyQAM() {}
  //! Class Constructor
  MyQAM(int M) { set_M(M); }
  //! Destructor
  virtual ~MyQAM() { }
  //! Change the size of the signal constellation
  virtual void set_M(int M);

  //! Hard demodulation of bits
  virtual void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  virtual bvec demodulate_bits(const cvec& signal) const;

protected:
  //! The square-root of M
  int L;
  //! Scaling factor of square QAM constellation (sqrt((M-1)*2/3))
  double scaling_factor;
};

class ITPP_EXPORT Paper19QAM : public MyQAM
{
public:
  //! Class Constructor
  Paper19QAM(): MyQAM() {set_M();}
  //! Destructor
  virtual ~Paper19QAM() {}

  void set_M();
  //! Hard demodulation of bits
  void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  bvec demodulate_bits(const cvec& signal) const;
};

class ITPP_EXPORT TriSq19QAM : public MyQAM
{
public:
  //! Class Constructor
  TriSq19QAM(): MyQAM() {set_M();}
  //! Destructor
  virtual ~TriSq19QAM() {}

  void set_M();
  //! Hard demodulation of bits
  void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  bvec demodulate_bits(const cvec& signal) const;
};
