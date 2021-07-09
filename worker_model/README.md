# worker_model
## Last update: 8 July 2021

Network model to explore the impact of working patterns and other work-related interventions on the transmission of SARS-COV-2 in a worker population.

## Associated papers
For details on the epidemiological model, network generation and parameterisation methods used, see:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009058

## Code usage
A full pipeline of the model, showing all functions, their locations and when they are called, is provided in 'worker_model_pipeline.pdf'. All functions also have associated docstrings for more information: type '?function_name' into the REPL to see docstring (after loading files into environment).

### Directory structure
#### include_files_network_model
Contains the bulk of the code, split into separate files according to purpose (defined at the top of each file).

#### cluster_shell_files
Files required to run code remotely on cluster machines

### Running the model
The main script for running the model is contained in 'worker_model.jl'. This can be run locally from the REPL, or from the command line. See script for details.

### Inputs
In order to run the model on the command line, there are 7 user-defined inputs (arguments) that are required:
1. JOB_ID
2. RNGSEED_BASE - along with JOB_ID, used to initialise random number generator
3. COUNTFINAL - number of replicates to run for each scenario combination
4. ENDTIME - length of each replicate (days)
5. CMAX - number of nodes (workers) in the network
6. WORKERTYPES - number of unique work sectors
7. RUNSETS - scenarios to run ([configurations, interventions])

If run locally from the REPL, the default values set in 'worker_model.jl' are used.

### Parameters and defaults
All parameters used in the model are defined in the 'parametertypes.jl' file. They are split into different structures depending on what they relate to.

Default values for most parameters are defined in the 'parametertypes.jl' file. Only parameters within the NetworkParameters and WorkplaceGenerationParameters structures have their defaults set outside this file, within the 'find_network_parameters()' function.

The network has been parameterised for the default 41 work sectors. It can be run for other numbers of sectors, however the defaults used (e.g. transmission risks and contact degree distributions in each sector) are arbitrarily set and should be revised before relying on any results.

### Configurations
A 'configuration' is used to set non-default parameter values for the generation of the network and simulation of disease spread. These parameters apply from, and are the same at, the start of every replicate performed with that configuration. These parameters can be changed during a replicate by using an 'intervention'.

A configuration can change any parameter that is contained within the ConfigurationVariables parameter structure. This structure contains 'copies' of parameters that are defined within other structures. When also defined in the ConfigurationVariables structure, its counterpart is automatically changed.

Configurations are defined within the 'load_configurations()' function, and given a name (string type). When running the model, the name of the desired configuration is given as an input. The "default" configuration does not change any parameters from their default values. A single configuration name may contain multiple different configurations (e.g. "adherence_svty" contains multiple configurations, each with a different value for the adherence parameter), and multiple configuration names can be given as input. When the model is run, it performs COUNTFINAL replicates of each configuration within all configuration names given.

New configurations can be defined by adding them to the 'load_configurations()' function. A configuration must have a defined 'n_configs' value that is greater than 0, defining the number of different configurations that are contained under that name.  All variables defined within the configuration (in order to change a parameter from its default value) must either have length = 1 (same across all n_configs configurations) or length = n_configs (changes within each configuration).

### Interventions
An 'intervention' is used to alter the parameters of the network and/or disease spread at a specified time during a replicate. This allows control measures to implemented, changed and removed at different times during the outbreak. We call the set of interventions used in a single simulation an 'intervention set'.

An intervention can change any parameter that is contained within the InterventionVariables parameter structure (in a similar way to configurations). These changes are reversed between replicates, so that the parameters given by the configuration apply at the start of every replicate.

Interventions are defined within the 'load_interventions()' function and given a name (string type). When running the model, each configuration name is paired with an intervention name in the inputted arguments (can be "none" if no interventions are required). As with configurations, a single intervention name can contain multiple intervention sets. Every intervention set is run for COUNTFINAL replicates on every configuration. For example, if n_configs = 10, n_intervention_sets = 10 and COUNTFINAL = 10, a total of 1000 replicates will be run.

New interventions can be defined by adding them to the 'load_interventions()' function. Every intervention set must be an array of InterventionVariables objects, one object for each intervention in the set. Each intervention must have a defined 'start_time'. Optionally, interventions can be reversed at a specified time by defining a 'reset_time' within the object.

Care should be taken when adding new parameters to the InterventionVariables parameter structure. Just because a parameter (or 'setting') is changed does not guarantee it will carry through to the network 'component'. See below.

### Settings and components
We sometimes refer to parameters as 'settings' and the parts of the network they influence 'components'.

Examples:
1. the setting 'network_generation_method_dynamic_social' is used to generate the component 'workday_social_contacts_by_day'
2. the settings 'sameday', 'ton' and 'toff' (defining work schedules) are used to define the component 'atwork'

Therefore, if these settings are changed, the related component must be regenerated. Regeneration functions have been written (in 'intervention_fns.jl'), and are automatically called (in 'affect_intervention!()'), for all settings currently in the InterventionVariables structure (see pipeline for a full list). New additions may require their own regeneration functions to be added to the 'intervention_fns.jl' file and called within the 'affect_intervention!()' function.

### Stochasticity
Some settings / components within the model are randomly sampled. These can either be sampled between every configuration or every replicate:

**Every replicate**
- atwork (regular weekly work schedule of nodes)
- adherence (including engagement with contact-tracing) and delay to adherence
- infection related times (e.g. time spent in latent phase, delay to symptom onset etc)
- returned to work status of each node (vs working from home)
- probability and relative infectiousness of asymptomatics
- initial disease states of nodes
- individual transmission risks in all settings
- occurrence of false negatives
- ability to recall infector (for backwards contact-tracing)

**Every configuration**
- All contact layers (even if changed by intervention, changes will always be the same and are reverted back to the same state between replicates)
- Workplace sizes and allocation of workers

### Plotting
The results can be plotted using the 'plot_worker_model.jl' script, following the same syntax used to run the model (note that the model must first be run and results saved - save and load file locations will need to be changed to your own local directory).

This script is currently set to plot the total weekly infections over time, for all configurations, intervention sets and replicates. This can be changed as necessary.

## Future improvements
- Triggered interventions (currently only simple examples available)
- Include children, schools and age structure
- More accurately parameterised healthcare and education sectors
- Allow different working schedules for different nodes (e.g. part-time and full-time)
- Calibration to data from different stages of the outbreak
- Profiling for speed and memory improvements
