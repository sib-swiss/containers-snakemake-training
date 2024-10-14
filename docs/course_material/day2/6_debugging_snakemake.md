## Additional advanced concepts

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/6_additional_concepts.pdf){: .md-button }

## Designing a workflow

There are many ways to design a new workflow, but these few pieces of advice will be useful in most cases:

* Start with a pen and paper: try to find out how many rules you will need and how they depend on each other. In other terms, start by sketching the DAG of your workflow!
    * Remember that Snakemake has a bottom-up approach (it goes from the final outputs to the first input), so it may be easier for you to work in that order as well and write your last rule first
    * Determine which rules (if any) aggregate or split inputs and create input functions accordingly (the topic is tackled in [this](5_reproducibility_snakemake.md) series of exercises)
* Make sure your input and output directives are right before worrying about anything else, especially the shell sections. There is no point in executing commands with the wrong inputs/outputs!
    * Remember that Snakemake builds the DAG before running the shell commands, so you can use the `-n/--dry-run/--dryrun` parameter to test the workflow before running it. You can even do that without writing all the shell commands!
* List any parameters or settings that might need to be adjusted later
* Choose meaningful and easy-to-understand names for your rules, inputs, outputs, parameters, `wildcards`... to make your Snakefile as readable as possible. This is true for every script, piece of code, variable... and Snakemake is no exception! Have a look at [The Zen of Python](https://peps.python.org/pep-0020/) for more information

## Debugging a workflow

It is very likely you will see bugs and errors the first time you try to run a new Snakefile: don’t be discouraged, this is normal!

### Order of operations in Snakemake

The topic was approached when we discussed [DAGs](2_introduction_snakemake.md#chaining-rules), but to efficiently debug a workflow, it is worth taking a deeper look at what Snakemake does when you execute the command `snakemake -c 1 <target>`. There are 3 main phases:

1. Prepare to run:
    1. Read all the rule definitions from the Snakefile
1. Resolve the DAG (happens when Snakemake says 'Building DAG of jobs'):
    1. Check what output(s) are required
    1. Look for matching input(s) by looking at the outputs of all the rules
    1. Fill in the `wildcards` to determine the exact input(s) of the matching rule(s)
    1. Check whether this(these) input(s) is(are) available; if not, repeat Step 2 until everything is resolved
1. Run:
    1. If needed, create the output(s) folder path
    1. If needed, remove the outdated output(s)
    1. Fill in the placeholders and run the shell commands
    1. Check that the commands ran without errors and produced the expected output(s)

### Debugging advice

Sometimes, Snakemake will give you a precise error report, but other times... less so. Try to identify which phase of execution failed (see previous paragraph on [order of operations](6_debugging_snakemake.md#order-of-operations-in-snakemake)) and double-check the most common error causes for that phase:

1. Parsing phase failures (phase 1):
    * Syntax errors, among which (but not limited to):
        * Missing commas/colons/semicolons
        * Unbalanced quotes/brackets/parenthesis
        * Wrong indentation
        * _These errors can be easily solved using a text editor with Python/Snakemake text colouring_
    * Failure to evaluate expressions
        * Problems in functions (`expand()`, input functions...) in input/output directives
        * Python logic added outside of rules
    * Other problems with rule definition
        * Invalid rule names/directives
        * Invalid wildcard names
        * Mismatched `wildcards`
1. DAG building failures (phase 2, before Snakemake tries to run any job):
    * Failure to determine the target
    * Ambiguous rules making the same output(s)
    * On the contrary, no rule making the required output(s)
    * Circular dependency (violating the **'Acyclic'** property of a D**A**G).
    * Write-protected output(s)
1. DAG running failures (phase 3, `--dry-run` works and builds the DAG, but the jobs execution fails):
    * Shell command returning non-zero status
    * Missing output file(s) after the commands have run
    * Reference to a `$shell_variable` before it was set
    * Use of a wrong/unknown placeholder inside `{ }`
    * _When a job fails, Snakemake reports an error, deletes all output file(s) for that job (because of potential corruption), and stops_
