## Learning outcomes

**After having completed this chapter you will be able to:**

* Understand the structure of a Snakemake workflow
* Write rules and Snakefiles to produce the desired outputs
* Chain rules together
* Run a Snakemake workflow

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/2_introduction_snakemake.pdf){: .md-button }

## Structuring a workflow

It is advised to implement your code in a directory called `workflow` (you will learn more about workflow structure in the [next series of exercises](3_generalising_snakemake.md#downloading-data-and-setting-up-folder-structure)). Filenames and locations are up to you, but we recommend that you at least group all workflow outputs in a `results` folder.

## Exercises

This series of exercises will bear no biological meaning, on purpose: it is designed to explain the fundamentals of Snakemake.

### Creating a basic rule

A **rule** is the smallest block of code with which you can build a **workflow**. It is a **set of instructions** to create one or more **output(s)** from zero or more **input(s)**. When a rule is **executed** (in other words, applied to specific input/output file(s)), it is called a **job**. The definition of a rule always starts with the **keyword** `rule`. Similarly to Python classes and their attributes, rules have **directives**, which contain information about their **properties**.

To create the simplest rule possible, you need at least two _directives_:

* `output`: path of the output file
* `shell`: shell commands that will create the output when they are executed

Other directives will be explained throughout the course.

**Exercise:** The following example shows the minimal syntax to implement a rule. What do you think it does? Does it create a file? If so, how is it called?

```python linenums="1"
rule hello_world:
    output:
        'results/hello.txt'
    shell:
        'echo "Hello world!" > results/hello.txt'
```

??? success "Answer"
    This rule uses the `echo` shell command to print `Hello world!` in an **output file** called `hello.txt`, located in the `results` folder.

Rules are defined and written in a file called **Snakefile** (note the capital `S` and the absence of extension in the filename). This file should be located at the workflow root directory (here, `workflow/Snakefile`).

### Executing a workflow with a specific output

It is now time to execute your first workflow! To do this, you need to tell Snakemake what is your **target**, _i.e._ what is the **specific output** that you want to generate. A target can be any output from any rule in the workflow.

**Exercise:** Create a Snakefile and copy the previous rule in it. Then, execute the workflow with `snakemake -c 1 <target>`. What value should you use for `<target>`? Once Snakemake execution is finished, can you locate the output file?

??? info "What does `-c/--cores` do?"
    The `-c/--cores N` parameter controls the maximum number of CPU cores used in parallel. If N is omitted or 'all', Snakemake will use all available CPU cores, which is useful but can also be dangerous on a cluster or a local machine. In case of cluster/cloud execution, this argument sets the maximum number of cores requested from the cluster or cloud scheduler.

??? warning "Code indentation in Snakemake"
    As Snakemake is built on top of Python, proper **code indentation is crucial**. Wrong indentation often results in cryptic errors. We recommend using **indents of 4 spaces**, but here are two rules that should be followed at all times:

    * Do not mix space and tab indents in a file
    * Always use the same indent length

??? success "Answer"
    * The target value is the file you want to generate, here `results/hello.txt`. The command to execute the workflow is:
    ```sh
    snakemake -c 1 results/hello.txt
    ```
    * The output is located in the `results` folder. You can check the folder content with `ls -alh results/`
    * You can check the output content with `cat results/hello.txt`

During the workflow execution, Snakemake automatically created the **missing folder** of the output path, `results/`. If several nested folders are missing (for example, `test1/test2/test3/hello.txt`), Snakemake will create the **entire folder structure** (`test1/test2/test3/`).

**Exercise:** Re-run the exact same command. What happens?

??? success "Answer"
    Nothing! You get a message saying that Snakemake did not run anything:

    ```python
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
    ```
    This is normal, because the desired output is already present and accounted for!

??? info "Snakemake re-run policy"
    By default, Snakemake runs a job if:

    * A target file explicitly requested in the `snakemake` command is missing or an intermediate file is missing and is required to produce a target file
    * It detects input files that have been modified more recently than output files, based on their modification dates. In this case, Snakemake will re-generate existing outputs
    * Code (including `params` directive, see [here](4_optimising_snakemake.md#non-file-parameters) for more information) has changed since last workflow execution
    * Computing environment has changed since last workflow execution

    Snakemake re-runs can be forced:

    * For a specific rule using the `-R/--forcerun` parameter: `snakemake -c 1 -R <rule_name>`
    * For a specific target using the `-f/--force` parameter: `snakemake -c 1 -f <target>`
    * For **all workflow outputs** using the `-F/--forceall` parameter: `snakemake -c 1 -F`

    In practice, Snakemake re-run policy can be altered, but we will not cover this topic in the course (see [--rerun-triggers parameter](https://snakemake.readthedocs.io/en/v8.20.5/executing/cli.html) in Snakemake CLI help and [this git issue](https://github.com/snakemake/snakemake/issues/1694) for more information).

In the previous rule, the values of the two directives are **strings**. In the `shell` directive (other types of values will be seen later in the course), long strings (which includes software commands) can be written on multiple lines for clarity by encasing each line in quotes:

```python linenums="1"
rule long_message:
    output:
        'results/long_message.txt'
    shell:
        'echo "I want to print a very very very very very '
        'very very very long string in my output" > results/long_message.txt'
```

Here, Snakemake will concatenate the two lines (_i.e._ paste the two lines together) and execute the resulting command:
```sh
echo "I want to print a very very very very very very very very very very long string in my output" > results/long_message.txt
```

### Understanding the input directive

Another directive used by most rules is `input`. It usually indicates a path to a file required by the rule to create the output. In the following example, we wrote a rule that uses the file `results/hello.txt` as an input, and copies its content to `results/copied_file.txt`:

```python linenums="1"
rule copy_file:
    input:
        'results/hello.txt'
    output:
        'results/copied_file.txt'
    shell:
        'cp results/hello.txt results/copied_file.txt'
```

You will use the `input` directive in the next exercises.

### Creating a workflow with several rules

As you may have guessed from the previous rule, the `input` and `output` directives allow us to create links (also called **dependencies**) between rules and files. Here, the `input` of rule `copy_file` requires the `output` of rule `hello_world`. In other terms, this is... **a workflow**! Let's build one with two rules and run it!

#### Rule order matters!

**Exercise:** Add the rule `copy_file` to your Snakefile, **after** rule `hello_world`. Then, run the workflow **without specifying an output** with `snakemake -c 1`. What happens?

??? tip "Your Snakefile should look like this"
    ```python linenums="1"
    rule hello_world:
        output:
            'results/hello.txt'
        shell:
            'echo "Hello world!" > results/hello.txt'

    rule copy_file:
        input:
            'results/hello.txt'
        output:
            'results/copied_file.txt'
        shell:
            'cp results/hello.txt results/copied_file.txt'
    ```

??? success "Answer"
    Nothing! You get the same message as before, saying that Snakemake did not run anything:

    ```sh
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
    ```

    When you do not specify a target, the one selected by **default** is the **output of the first rule in the Snakefile**, here `results/hello.txt` of rule `hello_world`. While this behaviour may seem weird, it will prove very useful later! In this case, `results/hello.txt` already exists from your previous runs, so Snakemake doesn't recompute anything.

Let's try to better understand how rule dependencies work in Snakemake.

#### Chaining rules

The execution principle behind Snakemake is to create a **Directed Acyclic Graph (DAG)** that defines **dependencies between all inputs and outputs** of the workflow. Starting from jobs generating the final desired outputs, Snakemake checks whether required inputs exist. If they do not, it looks for a rule that can generate these inputs and so on until all dependencies are resolved. This is why Snakemake is said to have a 'bottom-up' approach: it starts from last outputs and go back to first inputs.

??? bug "`MissingInputException`"
    `MissingInputException` is a common error in Snakemake. It means that Snakemake couldn't find a way to generate targets during DAG computation because an input file is missing. This is a case of **broken dependency** between rules. This error is often caused by typos in input or output paths (for example, output of rule `hello_world` not matching input of rule `copy_file`), so make sure to double-check them!

**Exercise:** With this in mind, identify the target you need to use to trigger the execution of rule `copy_file`. Add the `-F` parameter to the `snakemake` command and execute the workflow. What do you see?

??? info "What do we use `-F/--forceall` here?"
    The `-F/--forceall` parameter forces the re-creation of **all workflow outputs**. It is used here to avoid manually removing files, but it should be used carefully, especially with large workflows which contains a lot of outputs.

??? success "Answer"
    * To trigger the execution of the second rule, you need to use `results/copied_file.txt` as target. The command is:
    ```sh
    snakemake -c 1 -F results/copied_file.txt
    ```
    * You should now see Snakemake execute two rules and produce both targets/outputs: to generate output `results/copied_file.txt`, Snakemake requires input `results/hello.txt`. Before the workflow is executed, this file does not exist, therefore, Snakemake looks for a rule that generates `results/hello.txt`, here rule `hello_world`. The process is then repeated for `hello_world`. In this case, the rule does not require any input, so all dependencies are resolved, and Snakemake can generate the DAG

While it is possible to pass a space-separated list of targets in a Snakemake command, writing all the intermediary outputs does not look like a good idea: it is very time-consuming, error-prone... and annoying! Imagine what would happen with a workflow generating hundreds of files?! Using rule dependencies effectively solve this problem: you only need to ask Snakemake for the final outputs, and it will create the necessary intermediary outputs by itself!

#### Rule dependencies can be easier to write

Creating rule dependencies using long file paths can be cumbersome, especially when you are dealing with a large number of files/rules. But there is a dedicated Snakemake syntax that makes this process easier to set-up: it is possible (and recommended!) to refer to the output of a rule in another rule with the following syntax: `rules.<rule_name>.output`. It has several advantages, among which:

* It limits the risk of error because you do not have to write filenames in several locations
* Changes in the output name are automatically propagated to rules that use it, which means that you only need to change the name once, in the rule that defines it
* It makes the code much clearer and easier to understand: with this syntax, you instantly know the object type (a `rule`), how/where it is created (`hello_world`), and what it is (an `output`)

??? warning "Rules must produce unique outputs"
    Because of rule dependency, it is **mandatory** that an output be generated by a **single rule**. Rules generating the same output are called **ambiguous**. When Snakemake encounters ambiguous rules, it is not able to decide -at least by itself- which rule to use to generate this output, so it stops the execution. In reality, there are solutions to deal with ambiguous rules, but they should be avoided as much as possible, so we will not cover them in this course. See the [official documentation](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/rules.html#handling-ambiguous-rules) for more information).

??? info "To quote or not to quote?"
    As opposed to strings, like `'results/hello.txt'`, quotes are not required around `rules.<rule_name>.output` statements, because they are Snakemake objects.

The following example implements this syntax for the two rules defined above:

```python linenums="1" hl_lines="9"
rule hello_world:
    output:
        'results/hello.txt'
    shell:
        'echo "Hello world!" > results/hello.txt'

rule copy_file:
    input:
        rules.hello_world.output  # Dependency syntax
    output:
        'results/copied_file.txt'
    shell:
        'cp results/hello.txt results/copied_file.txt'
```

Try to use this syntax as much as possible in the next series of exercises!
