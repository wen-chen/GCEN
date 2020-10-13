#include <getopt.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "util/func.hpp"


void data_norm_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "data_norm usage:\n";
  std::cout << "  data_norm -i input_file -o output_file -m normalization_method\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -m  --method <upqt or median or deseq> normalization method (default: upqt)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m upqt\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    data_norm_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  std::string method = "upqt";

  const char * const short_opts = "hvi:o:m:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "method", 1, NULL, 'm' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        data_norm_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'o':
        out_file_name = optarg;
        break;
      case 'i':
        in_file_name = optarg;
        break;
      case 'm':
        method = optarg;
        break;
      case '?':
        data_norm_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long( argc, argv, short_opts, long_opts, NULL );
  }

  // check options
  if (in_file_name.empty()) {
    std::cerr << "Error: -i/--input is required but not specified!\n";
    return -1;
  }

  if (out_file_name.empty()) {
    std::cerr << "Error: -o/--output is required but not specified!\n";
    return -1;
  }

  // load file
  std::vector <std::string> RowName;
  std::vector <std::vector <double> > MatrixRow;
  load(in_file_name, RowName, MatrixRow);

  // matrix inversion
  unsigned int RowNum = MatrixRow.size();
  unsigned int ColNum = MatrixRow[0].size();
  std::vector <std::vector <double> > MatrixCol(ColNum, std::vector <double>(RowNum));
  matrix_transpose(MatrixRow, MatrixCol);

  // normalization
  std::vector <std::vector <double> > MatrixScaled(RowNum, std::vector <double>(ColNum));
  if (method == "upqt") {
    std::vector <double> upqts;
    for (unsigned int i = 0; i < ColNum; ++i) {
      std::vector <double> Col;
      for (unsigned int j = 0; j < RowNum; ++j) {
        if (MatrixCol[i][j] > 0) {
          Col.push_back(MatrixCol[i][j]);
        }
      }
      std::sort(Col.begin(), Col.end());
      double upqt = sorted_vector_upperquartile(Col);
      upqts.push_back(upqt);
    }

    double upqts_mean = 0.0;
    for (unsigned int i = 0; i < upqts.size(); ++i) {
      upqts_mean = upqts_mean + upqts[i];
    }
    upqts_mean = upqts_mean / upqts.size();

    std::vector <double> ratios;
    for (unsigned int i = 0; i < ColNum; ++i) {
      ratios.push_back(upqts_mean / upqts[i]);
    }

    for (unsigned int i = 0; i < RowNum; ++i) {
      for (unsigned int j = 0; j < ColNum; ++j) {
        MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
      }
    }
  } else if (method == "median") {
    std::vector <double> medians;
    for (unsigned int i = 0; i < ColNum; ++i) {
      std::vector <double> Col;
      for (unsigned int j = 0; j < RowNum; ++j) {
        if (MatrixCol[i][j] > 0) {
          Col.push_back(MatrixCol[i][j]);
        }
      }
      std::sort(Col.begin(), Col.end());
      double upqt = sorted_vector_median(Col);
      medians.push_back(upqt);
    }

    double medians_mean = 0.0;
    for (unsigned int i = 0; i < medians.size(); ++i) {
      medians_mean = medians_mean + medians[i];
    }
    medians_mean = medians_mean / medians.size();

    std::vector <double> ratios;
    for (unsigned int i = 0; i < ColNum; ++i) {
      ratios.push_back(medians_mean / medians[i]);
    }
        
    for (unsigned int i = 0; i < RowNum; ++i) {
      for (unsigned int j = 0; j < ColNum; ++j) {
        MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
      }
    }
  } else if (method == "deseq") {
    //对每一行的数据求几何平均数，如果该行的几何平均数不为零（即该行所有值都不为零），则该行每个数据除以几何平均数，得到一个比例数的矩阵。
    //比例数矩阵和原数据矩阵，列数一样。
    std::vector <std::vector <double> > MatrixRatio;
    for (unsigned int i = 0; i < RowNum; ++i) {
      double geo_mean = geometric_mean(MatrixRow[i]);
      if (geo_mean > 0) {
        std::vector <double> item;
        for (unsigned int j = 0; j < ColNum; ++j) {
          item.push_back((MatrixRow[i][j] / geo_mean));
        }
        MatrixRatio.push_back(item);
      }
    }
    //转置比例数矩阵
    unsigned int MatrixRatio_RowNum = MatrixRatio.size();
    unsigned int MatrixRatio_ColNum = MatrixRatio[0].size();
    std::vector <std::vector <double> > MatrixRatio_Col(MatrixRatio_ColNum, std::vector <double>(MatrixRatio_RowNum));
    for (unsigned int i = 0; i < MatrixRatio_RowNum; ++i) {
      for (unsigned int j = 0; j < MatrixRatio_ColNum; ++j) {
          MatrixRatio_Col[j][i] =  MatrixRatio[i][j];
      }
    }
    //对比例数矩阵的每一列求中位数。
    std::vector <double> medians;
    for (unsigned int i = 0; i < MatrixRatio_ColNum; ++i) {
      std::vector <double> Col;
      for (unsigned int j = 0; j < MatrixRatio_RowNum; ++j) {
        if (MatrixRatio_Col[i][j] > 0) {
          Col.push_back(MatrixRatio_Col[i][j]);
        }
      }
      std::sort(Col.begin(), Col.end());
      double median = sorted_vector_median(Col);
      medians.push_back(median);
    }
    //原数据矩阵的每一个值除以该列的(几何平均数的中位数)。
    double medians_mean = 0.0;
    for (unsigned int i = 0; i < medians.size(); ++i) {
      medians_mean = medians_mean + medians[i];
    }
    medians_mean = medians_mean / medians.size();
        
    std::vector <double> ratios;
    for (unsigned int i = 0; i < ColNum; ++i) {
      ratios.push_back(medians_mean / medians[i]);
    }

    for (unsigned int i = 0; i < RowNum; ++i) {
      for (unsigned int j = 0; j < ColNum; ++j) {
        MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
      }
    }
  }

  // count and output 
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }
  
  std::vector <double> mean_vec;
  std::vector <double> sd_vec;

  for (unsigned int i = 0; i < RowNum; ++i) {
    out_file << RowName[i];
    mean_vec.push_back(mean(MatrixScaled[i]));
    sd_vec.push_back(deviation(MatrixScaled[i]));
    for (unsigned int j = 0; j < ColNum; ++j) {
      out_file << '\t' << MatrixScaled[i][j];
    }
    out_file << '\n';
  }

  double mean_of_mean = mean(mean_vec);
  double sd_of_mean = deviation(mean_vec);
  std::cout << "The average gene expression mean is " << mean_of_mean << " +- " << sd_of_mean << ".\n";

  std::sort(mean_vec.begin(), mean_vec.end());
  double per5_mean = mean_vec[int (mean_vec.size() * 0.05)];
  double per10_mean = mean_vec[int (mean_vec.size() * 0.10)];
  double per20_mean = mean_vec[int (mean_vec.size() * 0.20)];
  double per25_mean = mean_vec[int (mean_vec.size() * 0.25)];
  double per33_mean = mean_vec[int (mean_vec.size() * 0.33)];
  double per50_mean = mean_vec[int (mean_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_mean << "; ";
  std::cout << "10%: " << per10_mean << "; ";
  std::cout << "20%: " << per20_mean << "; ";
  std::cout << "25%: " << per25_mean << "; ";
  std::cout << "33%: " << per33_mean << "; ";
  std::cout << "50%: " << per50_mean << ").\n";

  double mean_of_sd = mean(sd_vec);
  double sd_of_sd = deviation(sd_vec);
  std::cout << "The average gene expression sd is " << mean_of_sd << " +- " << sd_of_sd << ".\n";

  std::sort(sd_vec.begin(), sd_vec.end());
  double per5_sd = sd_vec[int (sd_vec.size() * 0.05)];
  double per10_sd = sd_vec[int (sd_vec.size() * 0.10)];
  double per20_sd = sd_vec[int (sd_vec.size() * 0.20)];
  double per25_sd = sd_vec[int (sd_vec.size() * 0.25)];
  double per33_sd = sd_vec[int (sd_vec.size() * 0.33)];
  double per50_sd = sd_vec[int (sd_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_sd << "; ";
  std::cout << "10%: " << per10_sd << "; ";
  std::cout << "20%: " << per20_sd << "; ";
  std::cout << "25%: " << per25_sd << "; ";
  std::cout << "33%: " << per33_sd << "; ";
  std::cout << "50%: " << per50_sd << ").\n";

  return 0;
}