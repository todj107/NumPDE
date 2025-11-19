# NumPDE

## How to use Github:

### Setup:
1) Ensure that you have a ssh key (Ask GPT).
2) Clone repo from: https://github.com/todj107/NumPDE
3) Ask for permission (Tor is the only one that can give that).

### When you have done changes:
1) *git add .* // Add everything in your folder
2) *git commit -m "Some message about your changes"*
3) *git pull* // Need to pull the new version on git
4) *git push* // Push your changes to git

### Trubles:
1) Never change the same file as another person. Git conflicts sucks.
2) If the terminal asks about your name or email. Read what it says and do that.
3) If the terminal askes if you want *rebase* or *no rebase* choose *rebase*

### Changed the wrong file (marge conflicts):
1) See witch files that are commited:  
	*git show --name-only*
2) Undo the commit, but keep the changes to the files:  
	*git reset --mixed HEAD~1*
3) Restore the file with unwanted changes to its previous state:  
	*git restore <file>*
4) Recommit:  
	*git add .*  
	*git commit -m "Conflict solved"*
5) Test that it worked (should not display the unwanted file):  
	*git show --name-only* 
6) Push/pull your commit:  
	*git pull*  
	*git push*

## Download required packeges:
In your environment run:
*pip install -r requirements.txt*


