If you run the first or third `R` scripts for the first time, models and some other time-consuming outputs will be saved as `RData` files here to be reused as needed, saving you some time later (you can bypass this behaviour to re-run models manually if you want; see code for details).

A copy of all the `R_output` folder containing all the `RData` we personally generated and used in the manuscript is available as a `zip` file under the "release" section of this repo (GitHub only, starting with release v1: https://github.com/mdahirel/pushed-pulled-2020-dynamics/releases )

Feel free to download it and unzip its content to your own `R_output` folder to save you some computing time. File structure should be as follow for files to be found by the scripts:

┗ 📂R_output
  ┣ 📂supplementary
  ┃ ┣ 📜model_S2_computercounts.Rdata
  ┃ ┣ 📜model_S2_computercountsbymacro.Rdata
  ┃ ┣ 📜model_S3_initialdispersal.Rdata
  ┃ ┣ 📜model_S7_sizefront_IBM.Rdata
  ┃ ┣ 📜model_S8_front_IBM.Rdata
  ┃ ┣ 📜model_S8_genetics_IBM.Rdata
  ┃ ┣ 📜model_S9_popsize.Rdata
  ┃ ┗ 📜README.md
  ┣ 📜loo_genet_expe.Rdata
  ┣ 📜model1_front_expe.Rdata
  ┣ 📜model2_genet_expe.Rdata
  ┣ 📜model2bis_genet_expe.Rdata
  ┣ 📜model2ter_genet_expe.Rdata
  ┣ 📜model3_front_IBM.Rdata
  ┣ 📜model4_genetics_IBM.Rdata
  ┗ 📜README.md  (a.k.a. the file you are reading right now)
