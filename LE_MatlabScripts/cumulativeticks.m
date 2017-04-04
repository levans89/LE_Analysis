% cumulativeticks()

ytickset{1,3}=ytickset{1,2};

for t = 2:length(ytickset)
    ytickset{t,3} = ytickset{t-1,3}+ytickset{t,2}
end