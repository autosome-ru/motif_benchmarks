cat config.json | ./jq -r '.ppm[] | join("\t")'
