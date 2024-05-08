Collaboration
~~~~~~~~~~~~~

The version releases for *SeidarT* are organized such that the first number is the version, the second number is a feature release/improvement, and the 3rd is a development release for beta testing. For example, 2.0.0 is a stable version for v2. A new feature addition would go to 2.1.0, and during the feature development, 2.0.* will be reserved for debugging and testing. Wheels for PyPi are automatically built on pull/push requests to *main* using `GitHub's workflows <https://docs.github.com/en/actions/using-workflows>`_. The YAML file for doing so can be found in the hidden directory *.github/workflows*. Any additional workflows need to be placed in this folder. 

Code Contribution
^^^^^^^^^^^^^^^^^

The GitHub repo is organized into 3 permanent branches - *main*, *development*, and *hotfix* - and temporary feature branches for isolating fixes and development. 


.. graphviz::

    digraph BranchFlow {
        graph [ranksep=1.75, nodesep=1.0, bgcolor="transparent", rankdir="TB"];
        node [shape=rectangle, fontname=Helvetica, fontsize=12];
        
        main [style="rounded,filled", fillcolor="#00ed4b", root=true];
        develop [style="rounded,filled", fillcolor="#ffb885"];
        hotfix [style="rounded,filled", fillcolor="#ff8f8f"];
        feature [style="rounded,filled", fillcolor="#ff6be1"];

        main -> develop [dir="both", label="New release cycle", fontcolor="#dcdcdc", color="#1cb82e", penwidth=2.5];
        develop -> feature [label="New feature", fontcolor="#dcdcdc", color=orange, penwidth=2.5];
        feature -> develop [label="Feature complete", color="#1cb82e", fontcolor="#dcdcdc", penwidth=2.5];
        develop -> hotfix [label="Issue in production", color="#c92037", fontcolor="#dcdcdc", penwidth=2.5];
        hotfix -> main [label="Issue resolved", color="#1cb82e", fontcolor="#dcdcdc", penwidth=2.5];
        hotfix -> develop [label="Changes merged", color="#1cb82e", fontcolor="#dcdcdc", penwidth=2.5];
        
        {rank=same; main develop}
        {rank=same; hotfix feature}
    }
   
Features should be pulled to the Development branch where they will be tested for release. Naming of feature branches needs to be in the form *feature::<name_of_feature>*

Issue Tracking
^^^^^^^^^^^^^^

For reporting bugs or requesting features, use the GitHub Issues to help the maintainers track development and for user reference. **Do not email the maintainer(s)**. Please provide enough information and keywords, and use labels to help enable accurate searching. Prior to filing an issue, search the issues to avoid duplication. 