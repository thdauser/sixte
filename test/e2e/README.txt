** Performing e2e tests with Sixte **

In order to perform a new e2e test with Sixte, the following steps
have to be taken:

- A folder with the suffix ".e2e" has to be created.

- This folder needs to contain an executable script called
  "run_test". The script can be written in any (installed) script
  language.

- Importantly, the script needs to either exit with "0" (success) or
  "1" (failure) and not crash (we need to be able to run all e2e
  tests, even if one test fails).

- A logfile and error file will be automatically created in this
  folder

- Testing can be done by running ./run_e2e_test [dir] manually

- General data and the setup is and should be stored in the data/
  directory
