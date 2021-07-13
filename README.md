## GATK4-PBSPRO-CROMWELL for HKU 
This repo should help colleagues at HKU to run workflows written in [wdl](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using a workflow management system ([Cromwell](https://cromwell.readthedocs.io/en/stable/)). 

- [test.conf](/test.conf) is a configuration file that helps us to familiarize with cromwell's default environment variables and how the engine works in the background. Use [test.conf](/test.conf) along with the [hello_world](/hello_world) directory, this is a great first step into learning wdl+cromwell.

- [my.conf](/my.conf) is the configuration file of the engine (cromwell). You can run jobs with or without Docker via Singularity using [my.conf](/my.conf), it is up to you, and whether you specifiy you want docker in the wdl scripts.

- [db.conf](/db.conf) is the version of [my.conf](/my.conf) that uses a different database than cromwell's default in memory database. This allows us to use call-caching i.e. detect when a job has been run in the past so that it doesn't have to re-compute results. This is a great way to save time if your pipeline crashes for whatever reason after a couple of hours/days. The successful tasks will be stored in the database so they don't need to be re-computed, you just need to compute again what crashed and everything that should come after the crash.

The other folders are ready to implement workflows which have been tested in our HPC.


------------

### Notes on requested pipelines [copy-number-variation](/copy-number-variation) and [mutect2](/mutect2). 

- It's important to highlight here that there is a great way to generate inputs for both pipelines. Since the pipelines themselves are constructed in a way to provide the maximum flexibility among the desired inputs (values) a pipeline has, every given pipeline takes a respective `json` file as input which can be filled according to the user needs.
- [WOMtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) is an excellent jar file application that parses a wdl file and provides a template with all the desired inputs in `json` format, this is the way to go to start creating your first inputs.
- You can download the `WOMtool` jar file from this [link](https://github.com/broadinstitute/cromwell/releases/tag/59).
- Once downloaded we can generate input templates for the `mutect2.wdl` pipeline in the following way: `java -jar womtool-59.jar inputs mutect2.wdl` , this will output to `stdout` a template with inputs a user can fill in values. To get `json` file we can try this: `java -jar womtool-59.jar inputs mutect2.wdl > mutect2.json`
- This is an [example](/mutect2/tmp.json) of how a template input looks like.

------

### Running workflows in high-throughput.

- Once you get the grasp on how pipelines are run in the Broad Institute, we can start executing them in a high-throughput manner. Cromwell is the execution engine that run wdl files, it requiores a backend configuration to run in our cluster (xomics). The configuration file used to get cromwell working is [my.conf](/my.conf). 
- Furthemore, cromwell can run in two different [modes](https://cromwell.readthedocs.io/en/stable/Modes/): run and server:
    - Run mode is a good way to get started with Cromwell and experiment quickly. Run mode launches a single workflow from the command line and exits `0` or `1` to indicate the result. 
    - Server mode is the mode you wish for most applications of Cromwell, suitable for production use (high-throughput). Server mode starts Cromwell as a web server that exposes REST endpoints. 

- In order to start a server you can use a create file like that and execute as a job in xomics:

    #!/bin/bash
    #PBS -N job_cromwell_server
    #PBS -l walltime=440:00:00
    #PBS -l select=1:ncpus=4:mem=20gb
    #PBS -q cgsd
    module load java
    java -Xmx15G -Dconfig.file=/groups/cgsd/$USER/gatk-workflows/my.conf -jar /groups/cgsd/$USER/cromwell-59.jar server

- After that you can issue the command `qstat -f` and look for the execution node in xomics on which your job is running. For example let's assume it runs in node `compute-0-3`. Once you know the node all you need is to create a `ssh` tunnnel with port forwarding between your local machine and the node where the server is running. You can accomlish that by executing the following command:

`ssh -N -f -L localhost:8000:compute-0-3:8000 $USER@xomics.cpos.hku.hk`

This command will forward a port (8000-default cromwell port) from you node to your local machine on port (8000) as well in this example. Once the port is opened all you need to do is to go to your browser and type http://localhost:8000/ and you will be presented with [Swagger](https://swagger.io/) interface used by cromwell to launch workflows in server mode.
Whenever you get to this point, you will notice that the Swagger interface exposes your local filesystem to launch workflows. For that I generally keep a copy of my github repository in my local machine since (`wdl` and `json`) files are very light and use them to launch them workflows in xomics. It's worth noting that all file paths present in the `json` file with the heavy (`bam`) files you may have stored in xomics are relative to the filesystem of xomics and Cromwell will use all filepaths relative to xomics, not to your local machine (which is great).

Finally if you wish to use Github and git as the distributed version control system for organizing and sharing code, I would advise to create your own repository on github or create a `branch` of the current one, so we can have different versions to compare against. I would request to not `commit` and `push` to the `main` branch directly if you choose to use the current repo.  

Happy coding ! :sunglasses: