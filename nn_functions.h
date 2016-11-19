#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstring>
#include <cmath>
#include <tuple>
#include <numeric>
#include <map>

#include "DTW.h"
#include "FastDTW.h"
#include "EuclideanDistance.h"

#define WINDOW_WIDTH 20

// Identity std::to_string
namespace std
{
    std::string to_string(std::string str)
    {
        return str;
    }
}
// Double quote a string
// alfa --> "alfa"
auto qs = [](std::string str)
{
    return "\"" + str + "\"";
};
// Create a key-value pair for json
// alfa, beta --> alfa : beta,\n
auto kv = [](std::string str, auto value, bool comma = true)
{
    cout << str << " : " << std::to_string(value) << (comma ? "," : "") << endl;
};
// Create an array or object for json
auto wrp = [](std::string pre, auto wrapped, std::string post, bool comma = false)
{
    cout << pre << endl;
    wrapped();
    cout << post << (comma ? "," : "") << endl;
};

using namespace fastdtw;

struct taggedTS {
    std::vector<double> ts_ret_data;
    std::vector<int> ts_abs_data;
    std::string ts_tag;
    int id;
    std::string UID;
};

double fastDTWdist (taggedTS query,
                    taggedTS candidate,
                    int use_time_domain,
                    int print_warp_path) {

    int query_end = query.ts_abs_data[query.ts_abs_data.size() - 1];
    int candidate_end = candidate.ts_abs_data[candidate.ts_abs_data.size() - 1];

    TimeSeries<double,1> tsI, tsJ;

    for (int i = 0; i<query.ts_ret_data.size(); ++i) {
        if (use_time_domain && query.ts_abs_data[i] > candidate_end) {
            break; // if the query goes on longer than the candidate, quit early
        }
        tsI.addLast(i, TimeSeriesPoint<double,1>(&(query.ts_ret_data)[i]));
    }

    for (int i = 0; i<candidate.ts_ret_data.size(); ++i) {
        if (use_time_domain && candidate.ts_abs_data[i] > query_end) {
            break; // if the candidate goes on longer than the query, quit early
        }
        tsJ.addLast(i, TimeSeriesPoint<double,1>(&(candidate.ts_ret_data)[i]));
    }

    if (tsI.size() == 0 || tsJ.size() == 0) {
        cout << "Timeseries of size 0 compared; exiting." << endl;
        abort();
    }

    TimeWarpInfo<double> info =
      FAST::getWarpInfoBetween(tsI,tsJ,WINDOW_WIDTH,EuclideanDistance());

    if (print_warp_path) {
        info.getPath()->print(std::cout);
    }

    return info.getDistance();
}

double fastDTWdist (taggedTS query,
                    taggedTS candidate,
                    int use_time_domain) {
    return fastDTWdist(query, candidate, use_time_domain, 0);
}

// TODO: perhaps find a better way to keep track of this.
int global_id = 0;

std::vector<taggedTS> load_TSfile(std::string fname, int verbose) {
    std::ifstream jobfile (fname.c_str(), std::ifstream::in);

    if (!jobfile) {
        cout << "No such file \"" << fname << "\" in folder." << endl;
        abort();
    }

    std::string tag_line;
    std::string ret_time_line;
    std::string abs_time_line;
    std::string tok;
    std::vector <taggedTS> tsbuffer;

    //get lines in groups of three - tag lines, return lines, abs time lines.
    while (std::getline(jobfile, tag_line)      &&
           std::getline(jobfile, ret_time_line) &&
           std::getline(jobfile, abs_time_line)) {

        taggedTS current_ts;

        // get the title and UID
        std::istringstream line_iss(tag_line);
        line_iss >> tok;
        current_ts.ts_tag = tok;
        line_iss >> tok;
        current_ts.UID = tok;

        // get the return times
        std::istringstream ret_iss(ret_time_line);
        while (ret_iss >> tok) {
            current_ts.ts_ret_data.push_back(std::atof(tok.c_str()));
        }

        // get the absolute times
        std::istringstream abs_iss(abs_time_line);
        while (abs_iss >> tok) {
            current_ts.ts_abs_data.push_back(std::atoi(tok.c_str()));
        }

        // set the taggedTS id.
        current_ts.id = global_id++;

        // push the completed taggedTS and increase the count.
        tsbuffer.push_back(current_ts);
    }

    if (verbose) {
        cout << "processed: " <<
            tsbuffer.size() <<
            " vectors in file " <<
            fname.c_str() << "\n";
    }

    return tsbuffer;
}

//compares query against dataset. returns 1 if
//predicted value is correct, 0 o/w.
void kNN_worker(taggedTS query,
                  std::vector<taggedTS> dataset,
                  std::vector<std::tuple<double, taggedTS>>& results,
                  int use_time_domain) {

	#pragma omp parallel for
    for (int i = 0; i < dataset.size(); ++i) {

        if (dataset[i].id == query.id) {
            cerr << "WARNING: query series is in reference space. Skipping it!\n";
            ++i;
            continue;
        }

        double this_result = fastDTWdist(query, dataset[i], use_time_domain);

		#pragma omp critical
        {
            results.emplace_back(this_result, dataset[i]);
        }
    }
}

// compares query against dataset.
void kNN_single(taggedTS query,
                std::vector<taggedTS> dataset,
                int use_time_domain) 
{
    // Vector of (distance, timeseries)
    std::vector<std::tuple<double, taggedTS>> results;
    // Run kNN, filling the above vector
    kNN_worker(query, dataset, results, use_time_domain);
    // Output the neighbours array
    wrp(qs("neighbours") + " : [", [&results]()
    {
        // Output formatter (outputs a single neighbor)
        auto outputter = [](bool last)
        {
            return [last](std::tuple<double, taggedTS> result)
            {
                double distance = std::get<0>(result);
                taggedTS neighbor = std::get<1>(result);
                wrp("{", [last, distance, neighbor]()
                {
                    kv(qs("distance"), distance);
                    kv(qs("tag"), qs(neighbor.ts_tag));
                    kv(qs("UID"), qs(neighbor.UID));
                    kv(qs("id"), neighbor.id, false);
                }, "}", last);
            };
        };

        // Output 0...n-1
        for_each(results.begin(), results.end()-1, outputter(false));
        // Output n
        outputter(true)(results.back());
    }, "]", true);
}

// compares query *list* against dataset.
void one_NN_many(std::vector<taggedTS> queryset, std::vector<taggedTS> dataset, int use_time_domain)
{
    // Output the outer array
    wrp("[ ", [&]()
    {
        // Output one element from the query-set
        auto outputter = [&](bool last)
        {
            return [&](taggedTS query)
            {
                wrp("{", [&]()
                {
                    // Output all neighbors of this query
                    kNN_single(query,
                               dataset,
                               use_time_domain);

                    // Output information on the query itself
                    wrp(qs("ground_truth") + " : {", [&]()
                    {
                        kv(qs("tag"), qs(query.ts_tag));
                        kv(qs("UID"), qs(query.UID));
                        kv(qs("id"), query.id, false);
                    }, "}");
                }, "}", last);
            };
        };

        // Output 0...n-1
        for_each(queryset.begin(), queryset.end()-1, outputter(false));
        // Output n
        outputter(true)(queryset.back());
    }, "]");
}
