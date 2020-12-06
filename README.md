## GATK4-PBSPRO-CROMWELL for HKU 
This repo should help colleagues at HKU to run workflows written in [wdl](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using a workflow management system (Cromwell). 
[my.conf](/my.conf) is the configuration file of the engine (cromwell). You can run jobs with or without Docker using [my.conf](/my.conf), it is up to you and wether you specifiy you want docker or not in the wdl scripts.
[db.conf](/db.conf) is the version of [my.conf](/my.conf) that uses a different database than cromwell's default in memory database. This allows us to use call-caching i.e. detect when a job has been run in the past so that it doesn't have to re-compute results.
