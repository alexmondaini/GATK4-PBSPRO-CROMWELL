## GATK4-PBSPRO-CROMWELL for HKU 
This repo should help colleagues at HKU to run workflows written in [wdl](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using a workflow management system (Cromwell). 

- [my.conf](/my.conf) is the configuration file of the engine (cromwell). You can run jobs with or without Docker using [my.conf](/my.conf), it is up to you and wether you specifiy you want docker or not in the wdl scripts.

- [db.conf](/db.conf) is the version of [my.conf](/my.conf) that uses a different database than cromwell's default in memory database. This allows us to use call-caching i.e. detect when a job has been run in the past so that it doesn't have to re-compute results.

- [test.conf](/test.conf) is a configuration file that helps us to familiarize with cromwell's default environment variables and how the engine works in the background. Use [test.conf](/test.conf) along with the [hello_world](/hello_world) directory, this is a great first step into learning wdl+cromwell.

The other folders are are ready to implement workflows which have been tested in our HPC. You don't need to modify anything there, they are to go !


Happy coding ! :sunglasses:
