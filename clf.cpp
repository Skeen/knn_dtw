#include "cmdline.h"
#include "nn_functions.h"

int main(int argc, char** argv) {
    struct gengetopt_args_info ai;

    if (cmdline_parser(argc, argv, &ai) != 0) {
        exit(1);
    }

    std::vector<taggedTS> query =
      load_TSfile(ai.query_filename_arg, ai.verbose_flag);
    std::vector<taggedTS> reference =
      load_TSfile(ai.reference_filename_arg, ai.verbose_flag);

    one_NN_many(query, reference, ai.use_time_domain_flag, ai.knn_arg);

    return 0;
}
