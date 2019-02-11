index.html: index.xml xml2html.rb citation.rb
	./xml2html.rb > $@
	xmllint --html --noout --valid $@

.PHONY: test
test: index.html
	xmllint --html --noout --valid $<

