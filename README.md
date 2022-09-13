##  GATK4 for Hong Kong University
This repo is designed to help HKU colleagues to run workflows written in [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) using [Cromwell](https://cromwell.readthedocs.io/en/stable/) as the workflow management system. Its is highly desirable that one has some familiarity with WDL and Cromwell before proceeding here.


Checking if your pipeline has the correct syntax, and creating inputs from WDL files:

------------
- [WOMtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) is an excellent jar file application that parses a wdl file and provides a template with all the desired inputs in `json` format, this is the way to go to start creating your first inputs.

- You can download the `WOMtool` jar file from this [link](https://github.com/broadinstitute/cromwell/releases/). Once downloaded we can generate input templates for the `mutect2.wdl` pipeline in the following way: `java -jar womtool-${version}.jar inputs mutect2.wdl` , this will output to `stdout` a template with inputs a user can fill in values. To get a `json` file we can do this: `java -jar womtool-59.jar inputs mutect2.wdl > mutect2.json`
------------


### Running workflows in high-throughput.

- Once you get the basics on how pipelines are run in WDL+Cromwell, we can start executing them in a high-throughput manner. Cromwell is the execution engine that run wdl files, it requires a backend configuration to run in our cluster (xomics). The configuration files to get Cromwell working in our cluster are:

- [test.conf](/test.conf) is a configuration file that helps us to familiarize with cromwell's default environment variables and how the engine works in the background. Use [test.conf](/test.conf) along with the [hello_world](/hello_world) directory, this is a great first step to learn and get started.

- [application.conf](/application.conf) is the production configuration file of the engine (cromwell), use this for production workflows.


- Furthemore, cromwell can be executed in two different [modes](https://cromwell.readthedocs.io/en/stable/Modes/). These are `run` and `server`:
    - `run` mode is a good way to get started with Cromwell and experiment quickly. Run mode launches a single workflow from the command line and exits with a `0` or `1` code to indicate the result. 
    - `server` mode is the mode you wish for most applications of Cromwell, suitable for production use (high-throughput). Server mode starts Cromwell as a web server that exposes REST endpoints. 

- In order to start `server` mode you can use a script such as the followig one:

```bash
#!/bin/bash
#PBS -N job_cromwell_server
#PBS -l walltime=440:00:00
#PBS -l select=1:ncpus=6:mem=60gb
#PBS -q cgsd
#PBS -k oe

# choose your cromwell version
VERSION=83
# choose your local port 
PORT=8000

cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -f -L ${PORT}:${HOSTNAME}:8000 ${USER}@xomics.cpos.hku.hk

   and point your web browser to http://localhost:8000
   
END


module load java/11.0.9
java -Xms10G -Xmx30G -Dconfig.file=/path/to/application.conf \
-jar /path/to/cromwell-${VERSION}.jar server
```

- Check the `stderr` returned from the script above when you send it as a job in the cluster, the `stderr` will have this name format `job_cromwell_server.eXXXX` where XXXX is your job number, copy the ssh line which has been expanded with the `hostname` and `port` and paste into your local machine terminal to create the tunnel.

- This command will forward a port (8000-default cromwell port) from your node to your local machine on port (8000) in this example. Once the port is opened, all you need to do is to go to your browser and type http://localhost:8000/ ,and you will be presented with the [Swagger](https://swagger.io/) interface used by cromwell to launch workflows in server mode.

- Whenever you get to this point, you will notice that the Swagger interface exposes your local filesystem to launch workflows. For that I keep a copy of this github repository in my local machine since (`wdl` and `json`) files are very light and use them to launch the workflows in xomics. It's worth noting that all file paths present in the `json` files are relative to the filesystem of xomics and Cromwell will use all filepaths relative to xomics, not to your local machine (which is great) and allows you to store the heavy stuff in the cluster and not in your local machine.

Finally if you wish to use Github and git as the distributed version control system for organizing and sharing code, I would advise to create your own repository on github or create a `branch` of the current one, so you can have different versions to compare against.

Happy coding ! :sunglasses:
