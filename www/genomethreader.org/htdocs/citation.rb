require 'rexml/document'

Citation = Struct.new("Citation",:author,:url,:title,:journal,:volume,:number,
                                 :pages,:year,:comment,:key)

def citation_iterator(filename)
  citation_list = Array.new()
  File.open(filename,"r") do |f|
   xml = REXML::Document.new(f)
    xml.elements.each do |citations|
      citations.elements.each do |citxml|
        citation = Citation.new(nil,nil,nil,nil,nil,nil,nil,nil,nil,nil)
        citxml.elements.each do |tag|
          case tag.name
            when "author"
              citation.author = tag.text
            when "url"
              citation.url = tag.text
            when "title"
              citation.title = tag.text
            when "journal"
              citation.journal = tag.text
            when "volume"
              citation.volume = tag.text
            when "number"
              citation.number = tag.text
            when "pages"
              citation.pages = tag.text
            when "year"
              citation.year = tag.text.to_i
            when "comment"
              citation.comment = tag.text
            when "key"
              citation.key = tag.text
          end
        end
        yield citation
      end
    end
    f.close_read
  end
end
