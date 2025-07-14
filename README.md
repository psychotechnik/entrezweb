# NCBI Entrez database Web interface, API and CLI client

## This projects makes use of the following libraries

REST API (django, django-ninja)
Web Frontend (htmx, hyperscript.js, flowbite)
CLI Client (typer, questionary, rich)
Background tasks (celery)

A request to download the sequence is initiated by the "download nucleotide" form in the web frontend and by the "download" command from the CLI. In both cases a celery task is created to process the request.

The sequences are downloaded as fasta files and saved to the filesystem. In the Django Nucleotide model they are represented as FileField instances.

The web frontend renders a table for a downloaded nucleotide. In order to do so it needs to parse the sequence into equal parts for each of the five columns in the table. The number of rows that are rendered in the table has to be limited since sequences can be very long. I split the string into parts using the textwrap.TextWrapper class. For the larger sequence, from part 2 of the assignment, textwrap takes too long resulting in the browser timing out. In order to avoid that a substring from the sequence needs to be split to prevent the timeout.

The search input field in the Web Frontend allows the user to submit a query using a regular expression. To parse the query I am using Python's re.finditer function. The matches are highlighted in the resulting table.

All requests for sequences from both the Web Frontend and the CLI application are made using a RESTful API powered by django-ninja. I chose that library because it requires less boilerplate in setting up and provides automatic API Docs.

The CLI client allows the user to download and list downloaded sequences. The client renders a table very similar to the Web Frontend table displaying 5 columns. The CLI Client displays a progress bar as the sequences are downloaded.

If I had more time I would implement a call to the esearch endpoint to allow the user to locate Entrez UIDs by species or some other metadata before they download the sequences.

Another improvement could improve performance could be caching search results or sequence rows in a key value store like redis or memcached.

Also, compressing the contents of files containing the sequences could be another improvement.

The API and the Web Frontend are available at
<https://entrez-web.kalinsky.me>

API Docs: <https://entrez-web.kalinsky.me/api/docs#/default/>
