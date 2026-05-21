## DQ streaming repository

- Download locally the AnalysisResults.root files per streaming production
  ```ruby
  python utils.py config.yml --download
  ```
- Merge locally the AnalysisResults.root files and get them per run (the pre-merging files can be removed)
  ```ruby
  python utils.py config.yml --merge
  ```
- Run the streaming_qa.C macro to produce the QA plots selecting the workflow (table-maker or table-reader)
  ```ruby
  root -l 'streaming_qa.C("table-maker")'
  ```
