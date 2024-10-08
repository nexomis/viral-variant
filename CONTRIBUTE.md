# Nextflow pipeline conventions at Nexomis

## 0. Workflow file structure
The file structure to follow for a nextflow project is as follows:

```
.
└── <worflow_name>/
    ├── conf/
    │   ├── ext.conf
    │   ├── publish.conf
    │   └── resources.conf           (*)
    ├── assets/                      (*)
    │   └── input_schema.json        (*)
    ├── modules/
    │   ├-- config/
    │   │   ├── process/
    │   │   │   └── labels.config
    │   │   └── profiles.config
    │   ├-- process/
    │   │   ├── <process_name_1>/
    │   │   │   └── main.tf
    │   │   └── <process_name_2>/
    │   │       └── main.tf
    │   └-- subworflow/
    │       └── <subworflow_name_1>/
    │           └── main.tf
    ├── main.nf
    ├── nextflow.config
    ├── nextflow_schema.json
    └── Readme.md

  '--'  : submodule
  '(*)' : optional
```

## 1. Help
The `help` will be managed by `plugin/nf-schema`:
  - Include in `main.nf`
  - Create the `nextflow_schema.json` file according to the parameter values in the `params{}` section in `nextflow.config`

## 2. Config
 will include the following config files:
Different configuration files can be embedded in `nextflow.config`, their inclusion order is fixed and is as follows (inverse to the order of priority): 
    - `modules/config/process/labels.config` (for managing resources assigned to processes/labels)
    - `modules/config/profiles.config`
    - `conf/publish.conf` (to manage the publication of processes and subworkflows)
    - `conf/ext.conf`
    - `modules/config/dag.config`
    - `modules/config/report.config`
    - `modules/config/timeline.config`
    - `conf/resources.conf`
**NOTE:** configuration parameters can be overridden by the `-c` option in Nextflow.

## 3. Modules

Except in very rare cases where no existing process or subworkflow calls are made, and where no creation of processes or subworkflows is relevant, the following conventions should be followed.

Since `include` works only at the repository level, modules will be linked in the `modules` folder (to avoid file duplication and prevent incompatibility issues during module updates).

```sh
mkdir modules
git submodule add https://github.com/nexomis/nf-subworkflows.git modules/subworkflows
git submodule add https://github.com/nexomis/nf-config.git modules/config
git submodule add https://github.com/nexomis/nf-process.git modules/process
```
**NOTE:** submodules are linked to a repository by their hash reference. For them to be updated, they need to be pulled/pushed again.

For each modularized element (subworkflow or subprocess), place it in a folder named after the element containing a `main.nf`, which will be sourced for calls.

### 4. Continue modularity: reusable, and easily maintainable

#### 4.1 Rules and Convention for Process
See [here](https://github.com/nexomis/nf-process/blob/main/README.md)

#### 4.2 Single Input Queue Channel
Although not a strict requirement, it is preferable that a (sub)workflow has only one queue channel as input, and that it be the first argument. If necessary, queue channels should be merged into tuples.

#### 4.3 Centralizing Publications
To easily adapt (without multiplying module versions) the publication of `processes` and `subworkflows` (not the case for `workflows`) specifically to each workflow, centralize publication operations in the `conf/publish.conf` file.  
**Note:** to target a specific call of a process or subworkflow they must be imported with an unique name. In all cases a process (or subworkflow) can't be call twices with the same name.

#### 4.4 Global Parameters only in main worflows
Do not directly call parameters (`params`) in processes or sub-workflows: processes and sub-workflows should not directly depend on global parameters. Use channels to pass parameter values in workflows: aggregate parameters into channels at the main workflow level, then pass these channels to processes and sub-workflows as inputs.
For steps that need to be skipped (usually `params.skip_xxxxx`), use branches on input channels (usefull to manage it by sample) or pass them through sub-workflow input parameters (e.g., a list of steps to skip).

#### 4.5 Typography
#### 4.5.1 Channel Instance
Use `lowerCamelCase` format.
#### 4.5.2 Groovy Variables (Other than Channels)
Use `snake_case` format.
#### 4.5.3 Groovy Functions
Use `lowerCamelCase` format.
#### 4.5.4 Processes or Sub-Workflows
Use `UPPER_SNAKE_CASE` format.
#### 4.5.5 Process Inputs and Outputs
Use `snake_case` format.
#### 4.5.6 Sub-Workflow Inputs and Outputs
Use `lowerCamelCase` format.

### 4.6 Handling Paths as Queue Channels
Unless in exceptional cases, paths in queue channels should be passed and retrieved as tuples: val(meta), path(file/dir). See details and exepection [here](https://github.com/nexomis/nf-process/blob/main/README.md#3-handling-paths-as-queue-channels)

## 5. Input files/directories by sample

When input files or directories need to be associated with specific samples in a run, it is recommended to use sample sheets. If the same input file or directory is referenced by multiple samples, Nextflow typically handles this with symbolic links rather than physically duplicating the files. However, this can still lead to redundant processing of the same file in different tasks, which can slow down the workflow.

To avoid this redundancy, it is advisable to use dedicated sample sheets for managing input files/directories without repeating processing. In these dedicated sample sheets, each input file or directory will be associated with a unique ID. This ID will then be used in the main sample sheet to reference the sample, instead of specifying the file/directory path directly.

This approach ensures that input files/directories are processed only once, even when linked to multiple samples.

Note: The sample sheets should be placed in the assets/ directory.

