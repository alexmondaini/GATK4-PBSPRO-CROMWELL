## GATK4-PBSPRO-CROMWELL for HKU 
This repo should help colleagues at HKU to run workflows written in [wdl](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using a workflow management system ([Cromwell](https://cromwell.readthedocs.io/en/stable/)). 

- [test.conf](/test.conf) is a configuration file that helps us to familiarize with cromwell's default environment variables and how the engine works in the background. Use [test.conf](/test.conf) along with the [hello_world](/hello_world) directory, this is a great first step into learning wdl+cromwell.

- [my.conf](/my.conf) is the configuration file of the engine (cromwell). You can run jobs with or without Docker via Singularity using [my.conf](/my.conf), it is up to you, and whether you specifiy you want docker in the wdl scripts.

- [db.conf](/db.conf) is the version of [my.conf](/my.conf) that uses a different database than cromwell's default in memory database. This allows us to use call-caching i.e. detect when a job has been run in the past so that it doesn't have to re-compute results. This is a great way to save time if your pipeline crashes for whatever reason after a couple of hours/days. The successful tasks will be stored in the database so they don't need to be re-computed, you just need to compute again what crashed and everything that should come after the crash.

The other folders are ready to implement workflows which have been tested in our HPC.


Happy coding ! :sunglasses:

------------

### Notes on requested pipelines [copy-number-variation](/copy-number-variation) and [mutect2](/mutect2). 

- It's important to highlight here that there is a great way to generate inputs for both pipelines. Since the pipelines themselves are constructed in a way to provide the maximum flexibility among the desired inputs (arguments) a pipeline has, every given pipeline takes a respective `json` file as input.
- [WOMtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) is an excellent jar file application that parses a wdl file and provides a template with all the desired inputs in `json` format.
- You can download the `WOMtool` jar file from this [link](https://github.com/broadinstitute/cromwell/releases/tag/59) 
- Once download we can generate input templates for the `mutect2.wdl` pipeline in the following way: `java -jar womtool-59.jar inputs mutect2.wdl` , this will output to `stdout` a template with inputs a user can fill in values. To get `json` file we can try this: `java -jar womtool-59.jar inputs mutect2.wdl > mutect2.json`
- This is an [example](/mutect2/tmp.json) of how a template input looks like.

------

### Running workflows in high-throughput.

Once you get the grasp on how pipelines are run in the Broad Institute, we can start executing them in a high-throughput manner.