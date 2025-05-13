function cOut = strsplit_mst(line, delimiter)

% Ein String wird anhand von Trennzeichen in Teilstrings zerlegt.
% Als Teilstring werden jeweils Zeichenfolgen erkannt, die durch ein
% Trennzeichen getrennt sind. Welche Zeichen als Trennzeichen gelten, wird
% durch den Parameter delimiter festgelegt.
%
% Syntax: cOut = strsplit(string, delimiter)
%
%   string ... String der zerlegt werden soll
%   delimiter ... ein/mehrere Trennzeichen
%                 z.B. {','} ... Beistrich als Trennzeichen
%                 z.B. {' ' , ';'} ... Leerzeichen und Semikolon
%                 z.B. {char(9)} ... Tabulator (ASCII-char-Code)
%                 z.B. {' ',char(9)} ... Tabulator und Leerzeichen
%
%   nrOfWords ... Anzahl der ermittelten Worte
%   cOut ... beinhaltet die Worte in Form von einzelnen Strings
%            z.B. cOut{1} = 'hello'
%                 cOut{2} = 'world'

cOut = [];
wordNr = 0;
charNr = 1;
while charNr <= length(line)
    
    if ~any(strcmp(line(charNr),delimiter))
        
        % new word detected
        wordNr = wordNr + 1;
        cOut{wordNr} = line(charNr);

        while charNr < length(line) % searching for end of word
            charNr = charNr+1;
            if any(strcmp(line(charNr),delimiter))
                break; % end of word
            else
                cOut{wordNr} = [cOut{wordNr} line(charNr)];
            end
        end
    end
    charNr = charNr + 1;
        
end