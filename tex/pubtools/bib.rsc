% no tabs as indentation, line size very big
print.align.key = 0
print.use.tab = off
print.line.length = 200

% sort order for fields
sort = on
%
check.rule { year "^[\"{]?[0-9][0-9][\"}]?$"   }
check.rule { year "^[\"{]?1[89][0-9][0-9][\"}]?$"   }
check.rule { year "^[\"{]?20[0-9][0-9][\"}]?$"   }
check.rule { year "" "\@ \$: Semantic error. Year has to be a suitable number" }

% Style Improvements.
%
% delete duplicate entries
check.double.delete = on

% delete empty fields
rewrite.rule {"^\" *\"$"}
rewrite.rule {"^{ *}$"}
rewrite.rule {"ˆ{}$"}
rewrite.rule {"ˆ\"\"$"}

% delete useless fields introduced by reference managers
delete.field {file}
delete.field {abstract}
delete.field {annote}
delete.field {keywords}

% correct page ranges
rewrite.rule {pages # "\([0-9]+\) *\(-\|---\) *\([0-9]+\)" = "\1--\3"}

