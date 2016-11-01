#define MAX_CLUSTERS 200
#define SIGMA_COEFF 2.0

#include "cmdline.h"
#include <search.h>
#include "nn_functions.h"

struct clusterDATA {
    float mean;
    float std_dev;
    std::string clus_tag;
};

int main(int argc, char** argv) {
    struct gengetopt_args_info ai;

    if (cmdline_parser(argc, argv, &ai) != 0) {
        exit(1);
    }

    int verbose = ai.verbose_flag;
    int use_time_domain = ai.use_time_domain_flag;
    char* que_FILENAME = ai.query_filename_arg;
    char* ref_FILENAME = ai.reference_filename_arg;

    std::vector<taggedTS> query =
      load_TSfile(que_FILENAME, verbose);
    std::vector<taggedTS> reference =
      load_TSfile(ref_FILENAME, verbose);

    if (verbose) {
        cout << "      ***** BEGIN CLUSTER MODEL FITTING*****" << endl << endl;
    }
    std::vector<taggedTS>::iterator reference_IT = reference.begin();
    std::vector<taggedTS>::iterator reference_IT_trail = reference.begin();
    std::vector<taggedTS>::iterator query_IT = query.begin();

    // Create a hashtable from tag -> cluster data
    hcreate(MAX_CLUSTERS);
    ENTRY tbl_insert_entry, tbl_query_entry;
    ENTRY* tbl_ret_entry;

    // Create a buffer of clusterDATA structs for the hashtable
    clusterDATA clusters[MAX_CLUSTERS];

    int vecs_in_cur_cluster = 1;
    int cluster_buf_cursor = 0;

    for (int i = 1; i < reference.size(); ++i) {
        reference_IT++;
        if (reference_IT_trail->ts_tag.compare(reference_IT->ts_tag) != 0) {
            cluster_std_dev(reference_IT_trail,
                            reference_IT_trail->ts_tag,
                            vecs_in_cur_cluster,
                            verbose,
                            clusters[cluster_buf_cursor].mean,
                            clusters[cluster_buf_cursor].std_dev,
                            use_time_domain);

            clusters[cluster_buf_cursor].clus_tag = reference_IT_trail->ts_tag;

            // load results into hashtable
            tbl_insert_entry.key = strdupa((reference_IT_trail->ts_tag).c_str());
            tbl_insert_entry.data = &(clusters[cluster_buf_cursor]);

            hsearch(tbl_insert_entry, ENTER);

            reference_IT_trail = reference_IT;
            vecs_in_cur_cluster = 1;
            ++cluster_buf_cursor;
        }
        else {
            ++vecs_in_cur_cluster;
        }
    }

    // there should be one last cluster
    cluster_std_dev(reference_IT_trail,
                    reference_IT_trail->ts_tag,
                    vecs_in_cur_cluster,
                    verbose,
                    clusters[cluster_buf_cursor].mean,
                    clusters[cluster_buf_cursor].std_dev,
                    use_time_domain);

    clusters[cluster_buf_cursor].clus_tag = reference_IT_trail->ts_tag;

    tbl_insert_entry.key = strdupa(reference_IT_trail->ts_tag.c_str());
    tbl_insert_entry.data = &(clusters[cluster_buf_cursor]);

    hsearch(tbl_insert_entry, ENTER);

    if (verbose) {
        cout << "      ***** END CLUSTER MODEL FITTING*****" << endl << endl;
        cout << "      ***** BEGIN CLASSIFICATION *****" << endl << endl;
    }

    // in case we want this; 1 greater than the cursor value.
    int total_clusters = cluster_buf_cursor + 1;
    reference_IT = reference_IT_trail = reference.begin();

    int true_knowns = 0;
    int true_unknowns = 0;
    int unknowns_reported = 0;
    int knowns_reported = 0;

    double best_dtw;
    taggedTS prediction;

    for (int i = 0; i < query.size(); ++i) {
        std::string ground_truth_tag = query_IT->ts_tag;

        one_NN_single(*(query_IT++),
                      reference,
                      best_dtw,
                      prediction,
                      use_time_domain);

        tbl_query_entry.key = strdupa((prediction.ts_tag).c_str());

        if (tbl_ret_entry = hsearch(tbl_query_entry, FIND)) {
            clusterDATA* this_cluster = (clusterDATA*) tbl_ret_entry->data;

            float mean_dist = std::abs(best_dtw - (this_cluster->mean));
            float confidence_radius = SIGMA_COEFF*(this_cluster->std_dev);

            if (mean_dist > confidence_radius) {
                ++unknowns_reported;
                if (verbose) {
                    cout << "prediction: UNKNOWN (nearest: " << prediction.ts_tag << ", UID: " << prediction.UID
                         << ") ground truth: " << ground_truth_tag << endl;
                }
                tbl_query_entry.key = strdupa(ground_truth_tag.c_str());

                // If the ground truth's tag is not a known cluster
                if (!hsearch(tbl_query_entry, FIND)) {
                    ++true_unknowns;
                }
            }
            else {
                ++knowns_reported;
                // this is the non-unknown branch
                if (verbose) {
                    cout << "prediction: " << prediction.ts_tag << " ground truth: "
                         << ground_truth_tag << " UID: " << prediction.UID << endl;
                }
                if (prediction.ts_tag.compare(ground_truth_tag) == 0) {
                    ++true_knowns;
                }
            }

            if (verbose) {
                cout << "    mean_dist was " << mean_dist
                     << " and class radius is " << confidence_radius << endl;
            }
        }
        else {
            cout << "ERROR: EXPECTED CLUSTERHASH NOT FOUND" << endl;
            abort();
        }

    }

    if (verbose) {
        cout << "      ***** END CLASSIFICATION*****" << endl << endl;

        cout << "      ***** BEGIN CLASSIFICATION STATS REPORT*****"
             << endl << endl;
    }

    int true_knowns_percent, true_unknowns_percent;

    if (knowns_reported > 0) {
        true_knowns_percent =
            int ((float) 100*true_knowns / knowns_reported);
    }
    else {
        true_knowns_percent = -1;
    }

    if (unknowns_reported > 0) {
        true_unknowns_percent =
            int ((float) 100*true_unknowns / unknowns_reported);
    }
    else {
        true_unknowns_percent = -1;
    }

    if (verbose) {
        cout << "Classifier queried on " << query.size() << " input vectors "
             << endl << "against " << reference.size() << " reference vectors. "
             << endl << endl;

        cout << total_clusters << " classes known in total." << endl << endl;

        cout << "Vectors more than " << SIGMA_COEFF
             << " standard deviations from their class mean"
             << endl << "were reported unknown." << endl << endl;

        cout << true_knowns << " correct predictions out of " << knowns_reported
             << ": " << true_knowns_percent
             << "\% accuracy." << endl;

        cout << true_unknowns << " true unknowns out of " << unknowns_reported
             << " reported" << ": " << true_unknowns_percent <<
             "\% accuracy." << endl;
    }
    else {
        cout << "{" << endl
             << "    " << "queries_file:'" << que_FILENAME << "'," << endl
             << "    " << "query_vector_count:" << query.size() << "," << endl
             << "    " << "reference_file:'" << ref_FILENAME << "'," << endl
             << "    " << "reference_vector_count:" << reference.size() << "," << endl
             << "    " << "unique_reference_classes:" << total_clusters << "," << endl
             << "    " << "std_dev_coefficient:" << SIGMA_COEFF << "," << endl
             << "    " << "true_knowns:" << true_knowns << "," << endl
             << "    " << "knowns_reported:" << knowns_reported << "," << endl
             << "    " << "true_knowns_percent:" << true_knowns_percent << "," << endl
             << "    " << "true_unknowns:" << true_unknowns << "," << endl
             << "    " << "unknowns_reported:" << unknowns_reported << "," << endl
             << "    " << "true_unknowns_percent:" << true_unknowns_percent << "," << endl
             << "}" << endl;
    }

    hdestroy();
    return 0;
}
