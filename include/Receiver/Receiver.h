#pragma once

// stl.
#include <memory>
#include <map>

// 3rd party.
#include <Eigen/Dense>

// forward decl.
class Options;

class Receiver {

  static long num;
  /** < Counter for number of receivers. */

  static std::vector<std::string> names;
  /** < Static vector holding all receiver names. */

  long mNum;
  /** < Number of this receiver. */

  double mPysLocX1, mPysLocX2, mPysLocX3;
  /** < Locations in the physical mesh */

  double mRefLocR, mRefLocS, mRefLocT;
  /** < Locations in element coordinates */

  std::string mName;
  /** < Receiver name */

 protected:

  std::map<std::string,std::vector<float>> store;
  /** < Store all fields in here */

 public:

  // See discussion I've had with myself in the Element header (mRec).
  virtual std::unique_ptr<Receiver> move() = 0;

  Receiver(std::unique_ptr<Options> const &options);
  virtual ~Receiver();
  static std::vector<std::shared_ptr<Receiver>> factory(std::unique_ptr<Options> const &options);

  inline void SetRefLocR (double val) { mRefLocR = val; }
  inline void SetRefLocS (double val) { mRefLocS = val; }
  inline void SetRefLocT (double val) { mRefLocT = val; }

  inline long Num() { return mNum; }
  inline std::string Name() { return mName; }
  inline double PysLocX1() { return mPysLocX1; }
  inline double PysLocX2() { return mPysLocX2; }
  inline double PysLocX3() { return mPysLocX3; }
  inline double RefLocR() { return mRefLocR; }
  inline double RefLocS() { return mRefLocS; }
  inline double RefLocT() { return mRefLocT; }

  void record(const double val, const std::string &field);

  virtual void write() = 0;

};
