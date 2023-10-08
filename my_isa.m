function answer = my_isa(x, classname)

answer = any(cellfun( @(x) strcmp(x, classname), cat(1, class(x), superclasses(x))));