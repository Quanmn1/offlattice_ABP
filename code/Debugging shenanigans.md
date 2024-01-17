# Jan 12 2024
The corner element of boxes suddenly turned into -1, while the corresponding particle didn't change box at all.

I noticed that it was right after a particle changed boxes in another region. So I looked carefully into the code changing the boxes, and discovered that I missed a conditional statement to check if an index is not -1 before changing that element of neighbors.

So I changed neighbors\[-2]. That turns out to be boxes\[4]\[4] because of the way memory in C works. My life is miserable. I'm done with C for today.