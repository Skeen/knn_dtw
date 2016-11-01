#include "cmdline.h"
#include "nn_functions.h"

int main(int argc, char** argv) {
    struct gengetopt_args_info ai;

    if (cmdline_parser(argc, argv, &ai) != 0) {
        exit(1);
    }

    /* yes, this is an ugly hack. just pass the same file
     * for both query and reference if both waves are in one file.
     */
    std::vector<taggedTS> query =
        load_TSfile(ai.query_filename_arg, ai.verbose_flag);
    std::vector<taggedTS> reference =
        load_TSfile(ai.reference_filename_arg, ai.verbose_flag);

    std::vector<taggedTS> joined;
    joined.reserve(query.size() + reference.size());
    joined.insert(joined.end(), query.begin(), query.end());
    joined.insert(joined.end(), reference.begin(), reference.end());

    taggedTS x,y;

    for (int i = 0; i < joined.size(); ++i) {
        if (joined[i].UID == ai.x_arg) {
            x = joined[i];
        }
        if (joined[i].UID == ai.y_arg) {
            y = joined[i];
        }
    }

    double dist = fastDTWdist(x, y, ai.use_time_domain_flag,
                              ai.print_warp_path_flag);

    cout << "\nCompared\n  " << x.ts_tag << ", " << x.UID << "\nto\n  " << y.ts_tag
         << ", " << y.UID << "\nand got DTW distance " << dist << ".\n";

    return 0;
}
