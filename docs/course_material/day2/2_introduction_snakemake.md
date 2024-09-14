## Learning outcomes

**After having completed this chapter you will be able to:**

* Understand the structure of a Snakemake workflow
* Write rules and Snakefiles to produce the desired outputs
* Chain rules together
* Run a Snakemake workflow

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/1_introduction_snakemake.pdf){: .md-button }

## Structuring a workflow

It is advised to implement your code in a directory called `workflow` (more information about this in the next series of exercises). You are free to choose the names and location of files for the different steps of your workflow, but, for now, we recommend that you at least group all outputs from the workflow in a `results` folder within the `workflow` directory.

## Exercises

This series of exercises will bear no biological meaning, on purpose: it is designed to explain the fundamentals of Snakemake.

### Creating a basic rule

Rules are the basic blocks of a Snakemake workflow. A **rule** is like a **recipe** indicating how to produce a specific **output** . The **execution** of a rule to create an output is called a **job**. A rule is defined in a Snakefile with the _keyword_ `rule` and contains _directives_ which indicate the rule's properties.

To create the simplest rule possible, we need at least two _directives_:

* `output`: path of the output file for this rule
* `shell`: shell commands to execute in order to generate the output

We will see other directives throughout the course.

**Exercise:** The following example shows the minimal syntax to implement a rule. What do you think it does? Does it create a file? If so, how is it called?

```python
rule first_step:
    output:
        'results/first_step.txt'
    shell:
        'echo "Hello world!" > results/first_step.txt'
```

??? success "Answer"
    This rule uses the `echo` shell command to print the line `Hello world!` in an **output file** called `first_step.txt`, located in the `results` folder.

Rules are defined and written in a file called **Snakefile** (note the capital `S` and the absence of extension in the filename). This file should be located at the root of the workflow directory (here, `workflow/Snakefile`).

### Executing a workflow with a specific output

It is now time to execute your first worklow! To do this, you need to tell Snakemake what is your **target**, *i.e.* what is the **specific output** that you want to generate. A target can be any output from any rule in the workflow.

**Exercise:** Create a Snakefile and copy the previous rule in it. Then, execute the workflow with `snakemake --cores 1 <target>`. What value should you use for `<target>`? Once Snakemake execution is finished, can you locate the output file?

??? warning "Code indentation in Snakemake"
    As Snakemake is built on top of Python, proper **code indentation is crucial**. Wrong indentation often results in cryptics errors. We recommend using **indents of 4 spaces**, but here are two rules that should be followed at all times:

    * Do not mix space and tab indents in a file
    * Always use the same indent length

??? success "Answer"
    * The target value is the file wou want to generate, here `results/first_step.txt`. The command to execute the workflow is `snakemake --cores 1 results/first_step.txt`
    * The output is located in the `results` folder. You can check the folder content with `ls -alh results/`
    * You can check the output content with `cat results/first_step.txt`

During the workflow execution, Snakemake automatically created the **missing folder** of the output path `results/`. If several nested folders are missing (for example, `test1/test2/test3/first_step.txt`), Snakemake will create the **entire folder structure** (`test1/test2/test3/`).

**Exercise:** Re-run the exact same command. What happens?

??? success "Answer"
    Nothing! We get a message saying that Snakemake did not run anything:

    ```python
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
    ```
    This is normal, because the desired output is already present and accounted for!

??? info "Snakemake re-run policy"
    By default, Snakemake only runs a job if:

    * A target file explicitly requested in the `snakemake` command is missing
    * An intermediate file is missing and is required to produce a target file
    * It detects input files that have been modified more recently than output files, based on their modification dates. In this case, Snakemake will re-generate the existing outputs

    Snakemake re-runs can be forced:

    * For a specific rule using the `-R` option: `snakemake --cores 1 -R <rule_name>`
    * For a specific target using the `-f` option: `snakemake --cores 1 -f <target>`
    * For **all workflow outputs** using the `-F` option: `snakemake --cores 1 -F`

    In practice, we can also alter Snakemake re-run policy, but we will not cover this topic in the course (see [--rerun-triggers option](https://snakemake.readthedocs.io/en/stable/executing/cli.html) in Snakemake's CLI help and [this git issue](https://github.com/snakemake/snakemake/issues/1694) for more information).

In the previous example, the values of the two rule directives are **strings**. For the `shell` directive (we will see other types of values later in the course), long strings (which includes software commands) can be written on multiple lines for clarity, simply using a set of quotes for each line:

```python
rule first_step:
    output:
        'results/first_step.txt'
    shell:
        'echo "I want to print a very very very very very very '
        'very very very very long string in my output" > results/first_step.txt'
```

Here, Snakemake will concatenate the two lines (_i.e._ paste the two lines together) and execute the resulting command (`echo "I want to print a very very very very very very very very very very long string in my output" > results/first_step.txt`).

### Understanding the input directive

Another directive used by most rules is the `input` directive. It usually indicates the path to a file that is required by the rule to generate the output. In the following example, we wrote a rule that uses the file `results/first_step.txt` as an input, and copies its content to `results/second_step.txt`:

```python
rule second_step:
    input:
        'results/first_step.txt'
    output:
        'results/second_step.txt'
    shell:
        'cp results/first_step.txt results/second_step.txt'
```

We will use the `input` directive in the next exercises.

### Creating a workflow with several rules

As you may have guessed from the previous rule, the `input` and `output` directives allow us to create links (or rather, **dependencies**) between rules and files. Here, the `input` of `rule second_step` requires the `output` of `rule first_step`. In other terms, this is... **a workflow**! Let's build one with two rules and run it!

#### Rule order matters!

**Exercise:** Add the rule `second_step` to your Snakefile, **after** the `first_step` rule. Then, run the workflow **without specifying an output** with `snakemake --cores 1`. What happens?

??? tip "Your Snakefile should look like this"
    ```python
    rule first_step:
        output:
            'results/first_step.txt'
        shell:
            'echo "Hello world!" > results/first_step.txt'

    rule second_step:
        input:
            'results/first_step.txt'
        output:
            'results/second_step.txt'
        shell:
            'cp results/first_step.txt results/second_step.txt'
    ```

??? success "Answer"
    Nothing! We get the same message as before, saying that Snakemake did not run anything:

    ```
    Building DAG of jobs...
    Nothing to be done (all requested files are present and up to date).
    ```

    When you do not specify a target, the one selected by **default** is the **output of the first rule in the Snakefile**, here `results/first_step.txt` of `rule first_step`. While this behaviour may seem weird, it will prove very useful later! In this case, `results/first_step.txt` already exists from your previous runs, so Snakemake doesn't recompute anything.

Now, let's try to understand how rule dependencies work in Snakemake.

#### Chaining rules

The core principle of Snakemake's execution is to compute a **Directed Acyclic Graph (DAG)** that summarizes **dependencies between all the inputs and outputs** required to generate the final desired outputs. Starting from the jobs generating the final outputs, Snakemake checks whether the required inputs exist. If they do not, it looks for a rule that generates these inputs. This process is repeated until all dependencies are resolved. This is why Snakemake is said to have a 'bottom-up' approach: it starts from the last outputs and go back to the first inputs.

??? bug "`MissingInputException`"
    The `MissingInputException` error is common in Snakemake. It means that Snakemake couldn't a find to generate the targets during the DAG computation because an input file is missing. This is a case of 'broken dependency' between rules. This error is often caused by typos in input or output paths (for example,  the output of rule `first_step` not matching the input of rule `second_step`), so make sure to double-check them!

**Exercise:** With this in mind, identify the target you need to use to trigger the execution of `rule second_step`. Add the `-F` option to the `snakemake` command and execute the workflow. What do you see?

??? info "What does `-F` do?"
    The `-F` option forces the re-creation of **all outputs** of the workflow. It is used here to avoid manually removing files, but it should be used carefully, especially with large workflows which contains a lot of outputs.

??? success "Answer"
    * To trigger the execution of the second rule, you need to use `results/second_step.txt` as target. The command is `snakemake --cores 1 -F results/second_step.txt`
    * You should now see Snakemake execute the two rules and produce both targets/outputs: to generate the output `results/second_step.txt`, Snakemake requires the input `results/first_step.txt`. Before the workflow is executed, this file does not exist, therefore, Snakemake looks for a rule that generates `results/first_step.txt`, here the rule `first_step`. The process is then repeated for `first_step`. In this case, the rule does not require any input, so all dependencies are resolved, and Snakemake can generate the DAG

While it is possible to pass a space-separated list of targets in a Snakemake command, writing all the intermediary outputs does not look like a good idea: it is very time-consuming, error-prone... and annoying! Imagine what would happen with a workflow generating hundreds of files?! Using rule dependencies effectively solve this problem: you only need to ask Snakemake for the final outputs, and it will create the necessary intermediary outputs by itself!

#### Rule dependencies can be easier to write

Creating rule dependencies using long filepaths can be cumbersome, especially when you are dealing with a large number of files/rules. But there is a dedicated Snakemake syntax that makes this process easier to set-up: it is possible (and recommended!) to refer to the output of a rule in another rule with the following syntax: `rules.<rule_name>.output`. It has several advantages, among which:

* It limits the risk of error because we do not have to write filenames in several locations
* A change in the output name will be automatically propagated to rules that depend on it, *i.e.* the name only has to be changed once
* It makes the code much clearer and easier to understand: with this syntax, we instantly know the object type (a `rule`), how/where it is created (`first_step`), and what it is (an `output`)

??? warning "Rules must produce unique outputs"
    Because of the rule dependency process, by default, an output can only be generated by a single rule. Otherwise, Snakemake cannot decide which rule to use to generate this output, and the rules are considered **ambiguous**. In practice, there are ways to deal with ambiguous rules, but they should be avoided as much as possible, so we will not cover them in this cours. See the [official documentation](https://snakemake.readthedocs.io/en/v8.20.3/snakefiles/rules.html#handling-ambiguous-rules) for more information).

??? info "To quote or not to quote"
    As opposed to strings, like `'results/first_step.txt'`, you don't need quotes around `rules.<rule_name>.output` statements, because they are Snakemake objects.

The following example implements this syntax for the two rules defined above:

```python
rule first_step:
    output:
        'results/first_step.txt'
    shell:
        'echo "Hello world!" > results/first_step.txt'

rule second_step:
    input:
        rules.first_step.output  # Dependency syntax
    output:
        'results/second_step.txt'
    shell:
        'cp results/first_step.txt results/second_step.txt'
```

Try to use this syntax as much as possible in the other series of exercises!
