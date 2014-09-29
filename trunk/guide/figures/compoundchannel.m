
open('compoundchannel.fig');
plot([0 0 1 1 1.25 1.25 2.25 2.25],...
     [0 -.1 -.1 -.5 -.5 -.1 -.1 0],'k-',...
     [0 2.25],[-.02 -.02],'k-');
axis image;
axis off;

print -deps2 compoundchannel