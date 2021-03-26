# Opinion-Dynamics
This code comes from my Master's level dissertation titled "The Mathematics of Opinion Formation", which was supervised by Dr Benjamin Goddard (University of Edinburgh) who helped me many times to solve problems and debug my code.

## What are opinion dynamics?
The study of opinion dynamics involves modelling how the opinions of a group of individuals change over time. Say we have a group of people who are discussing how they feel about marmite. The opinions will range from "I hate marmite" to "I love marmite", but these opinions could be changed by interacting with the rest of the group. The main assumption we make in this project is that people will only be influenced by opinions similar to their own, so someone who loves marmite is not likely to be influenced by someone who hates it. 

We also investigate the inclusion of groups known as radicals. A group of radicals have a fixed and extreme opinion and cannot be influenced by anyone else. In our marmite example, this is like introducing a small group of marmite haters into the main group. The marmite hating group then refuse to listen to the rest of the group, but will try to convince everyone else that they should hate marmite.

## What is included in this repo?
This repo contains code which implement three different opinion dynamics models, as well as code to compare the results of the SDE and PDE models. Within the SDE_model and PDE_model folders, you will also find testing folders which contain scripts used to assess the behaviour of our models.
