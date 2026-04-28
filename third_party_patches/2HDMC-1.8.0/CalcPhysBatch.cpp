#include "HBHS.h"
#include "THDM.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace {

bool build_point(double mh_in, double mH_in, double mA_in, double mHp_in,
                 double sba_in, double l6_in, double l7_in, double m12_2_in,
                 double tb_in, int yt_in, const string &output_file,
                 string &message) {
  THDM model;

  SM sm;
  model.set_SM(sm);

  bool pset = model.set_param_phys(mh_in, mH_in, mA_in, mHp_in, sba_in, l6_in,
                                   l7_in, m12_2_in, tb_in);
  if (!pset) {
    message = "invalid_parameters";
    return false;
  }

  model.set_yukawas_type(yt_in);

  const HBHSResult *hbhsres_ptr = nullptr;
#if defined HiggsBounds
  HBHS hbhs{};
  const HBHSResult hbhs_result = hbhs.check(model);
  hbhsres_ptr = &hbhs_result;
#endif

  model.write_LesHouches(output_file.c_str(), true, true, true, hbhsres_ptr);
  message = "ok";
  return true;
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: ./CalcPhysBatch input_table output_dir\n";
    cerr << "Each non-comment line must contain:\n";
    cerr << "row_id mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 m12^2 tan_beta yukawas_type\n";
    return -1;
  }

  ifstream input(argv[1]);
  if (!input) {
    cerr << "Could not open input table: " << argv[1] << "\n";
    return -1;
  }

  string output_dir = argv[2];
  string line;
  size_t processed = 0;
  size_t failed = 0;

  while (getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    istringstream stream(line);
    string row_id;
    double mh_in;
    double mH_in;
    double mA_in;
    double mHp_in;
    double sba_in;
    double l6_in;
    double l7_in;
    double m12_2_in;
    double tb_in;
    int yt_in;

    if (!(stream >> row_id >> mh_in >> mH_in >> mA_in >> mHp_in >> sba_in >> l6_in >> l7_in >> m12_2_in >> tb_in >> yt_in)) {
      cerr << "Malformed input line: " << line << "\n";
      failed += 1;
      continue;
    }

    string output_file = output_dir + "/" + row_id + ".lha";
    string message;
    bool ok = build_point(mh_in, mH_in, mA_in, mHp_in, sba_in, l6_in, l7_in,
                          m12_2_in, tb_in, yt_in, output_file, message);

    cout << row_id << '\t' << (ok ? "ok" : "failed") << '\t' << output_file
         << '\t' << message << "\n";
    processed += 1;
    if (!ok) {
      failed += 1;
    }
  }

  cerr << "Processed " << processed << " points, failed " << failed << "\n";
  return failed == 0 ? 0 : 1;
}