## GATK4-PBSPRO-CROMWELL for HKU 
This repo should help colleagues at HKU to run workflows written in [wdl](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using a workflow management system ([Cromwell](https://cromwell.readthedocs.io/en/stable/)). 

- [test.conf](/test.conf) is a configuration file that helps us to familiarize with cromwell's default environment variables and how the engine works in the background. Use [test.conf](/test.conf) along with the [hello_world](/hello_world) directory, this is a great first step into learning wdl+cromwell.

- [my.conf](/my.conf) is the configuration file of the engine (cromwell). You can run jobs with or without Docker via Singularity using [my.conf](/my.conf), it is up to you, and whether you specifiy you want docker in the wdl scripts.

- [db.conf](/db.conf) is the version of [my.conf](/my.conf) that uses a different database than cromwell's default in memory database. This allows us to use call-caching i.e. detect when a job has been run in the past so that it doesn't have to re-compute results. This is a great way to save time if your pipeline crashes for whatever reason after a couple of hours/days. The successful tasks will be stored in the database so they don't need to be re-computed, you just need to compute again what crashed and everything that should come after the crash.

The other folders are ready to implement workflows which have been tested in our HPC.


Happy coding ! :sunglasses:
