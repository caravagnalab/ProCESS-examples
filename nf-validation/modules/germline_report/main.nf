process GENERATE_REPORT {
    tag "generate_report"

    input:
    val caller_rds_map         // a map: caller → list of rds files
    path ground_truth_rds      // the ground truth .rds file

    output:
    path "dummy_report.txt", emit: report

    script:
    """
    echo "==== GERMLINE REPORT DUMMY ====" > dummy_report.txt
    echo "Ground truth file: ${ground_truth_rds}" >> dummy_report.txt
    echo "Caller-wise input RDS files:" >> dummy_report.txt

    ${caller_rds_map.collect { caller, rds_list -> 
        def rds_str = rds_list.collect { it.toString() }.join(' ')
        "echo '${caller}: ${rds_str}' >> dummy_report.txt"
    }.join('\n')}
    """
}

