library("devtools")

res <- revdep_check("inbreedR")
revdep_check_save_summary(res)
revdep_check_save_logs(res)
