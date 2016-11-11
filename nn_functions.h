#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstring>
#include <cmath>

#include "DTW.h"
#include "FastDTW.h"
#include "EuclideanDistance.h"

#define WINDOW_WIDTH 20

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
int one_NN_single(taggedTS query,
                  std::vector<taggedTS> dataset,
                  double& best_dtw,
                  taggedTS& prediction,
                  int use_time_domain) {

    best_dtw = 99999999;

	#pragma omp parallel for
    for (int i = 0; i < dataset.size(); ++i) {

        if (dataset[i].id == query.id) {
            cout << "WARNING: query series is in reference space. Skipping it!\n";
            ++i;
            continue;
        }

        double this_result =
          fastDTWdist(query, dataset[i], use_time_domain);

		#pragma omp critical
		{
        	if (this_result < best_dtw) 
			{
        		best_dtw = this_result;
        		prediction = dataset[i];
        	}
		}
    }

    return (query.ts_tag.compare(prediction.ts_tag) == 0);
}

//compares query *list* against dataset. returns
//number of correct classifications.
int one_NN_many(std::vector<taggedTS> queryset, std::vector<taggedTS> dataset, int use_time_domain) 
{
	//unsigned int concurrency = std::thread::hardware_concurrency();
	
   	int successes = 0;

	cout << "{" << endl;
	cout << " data : [" << endl;

	//#pragma omp parallel for reduction(+:successes)
    for (int i = 0; i < queryset.size(); ++i) 
	{
    	double best_dtw;
    	taggedTS prediction;
        successes += one_NN_single(queryset[i],
                                   dataset,
                                   best_dtw,
                                   prediction,
                                   use_time_domain);
		cout << "{" << endl;
		cout << "nearest_neighbor : '" << prediction.ts_tag << "'," << endl;
		cout << "UID : '" << prediction.UID << "'," << endl;
		cout << "distance : " << best_dtw << "," << endl;
		cout << "ground_truth : '" << queryset[i].ts_tag << "'," << endl;
		cout << "ground_truth_UID : '" << queryset[i].UID << "'" << endl;
		
		if(i == queryset.size() -1)
		{
			cout << "}"  << endl;
		}
		else
		{
			cout << "}," << endl;
		}
        //cout << "nearest neighbor: " << prediction.ts_tag << " (UID: " << prediction.UID
        //     << ") with distance " << best_dtw << ". ground: " << queryset[i].ts_tag
        //     << " (UID: " << queryset[i].UID << ")" << endl;
    }

    fflush(stdout);

	cout << "]," << endl;
	cout << "result : { " << endl;
	cout << "successes : " << successes << "," << endl;
	cout << "trials : " << queryset.size() << "," << endl;
	cout << "accuracy : " << ((float) 100*successes / queryset.size()) << endl;
	cout << "}" << endl;
	cout << "}" << endl;
    //cout << successes << " successes in " << queryset.size() << " trials: " << int ((float) 100*successes / queryset.size()) << "\% accuracy.\n";
    return successes;
	//return 0;
}

int cluster_std_dev(std::vector<taggedTS>::iterator cluster,
                    std::string clustag,
                    int vecs_in_cluster,
                    int verbose,
                    float& mean,
                    float& sigma,
                    int use_time_domain) {

    if (vecs_in_cluster < 2) {
        cout << "Clusters are size 2 or greater.\n";
        abort();
    }

    std::vector<taggedTS>::iterator cluster_IT = cluster;
    std::vector<double> distances;

    if (verbose) {
        cout << vecs_in_cluster << " vectors in cluster: " << clustag << endl;
    }

    for (int i = 0; i < vecs_in_cluster - 1; ++i) {
        std::vector<taggedTS>::iterator cluster_IT_i = cluster_IT;
        // we want to compare the timeseries at cluster_IT_i
        // with all the ones after it, but in the cluster.
        for (int j = i + 1; j < vecs_in_cluster; ++j) {
            distances.push_back(
                    fastDTWdist(*(++cluster_IT_i), *cluster_IT, use_time_domain));
        }
        ++cluster_IT;
    }

    // go back to head in case we want to use this again
    // for some reason later
    cluster_IT = cluster;

    float mean_accum = 0;
    for (int i = 0; i < distances.size(); ++i) {
        mean_accum += distances[i];
    }
    mean = mean_accum / distances.size();

    if (verbose) {
        cout << "Got mean: " << mean << endl;
    }

    float var_accum = 0;
    for (int i = 0; i < distances.size(); ++i) {
        var_accum += (mean - distances[i]) * (mean - distances[i]);
    }
    sigma = std::sqrt(var_accum / distances.size());

    if (verbose) {
        cout << "Got std.dev.: " << sigma << endl;
    }

    return 0;
}
