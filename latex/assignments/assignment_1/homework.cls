\NeedsTeXFormat{LaTeX2e}
\LoadClassWithOptions{article}
\ProvidesClass{homework}[Class file for homework assignments MATH301 2021]

% ----- Options ---------------------------------------------------------------
\newcommand\@opanon{0}
\DeclareOption{anonymous}{\renewcommand\@opanon{1}}
\newcommand\@opnewpage{0}
\DeclareOption{newpage}{\renewcommand\@opnewpage{1}}
\newcommand\@oplargemargins{0}
\DeclareOption{largemargins}{\renewcommand\@oplargemargins{1}}
\ProcessOptions

% ----- Packages --------------------------------------------------------------

% Better fonts with accents
\RequirePackage[T1]{fontenc}
%\RequirePackage{braket}

% Math symbols
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amsthm}
\RequirePackage{amssymb}
\RequirePackage{centernot}

% Nice lists
\RequirePackage{enumerate}
\RequirePackage{enumitem}

% Nice images, figures, and listings
\RequirePackage{graphicx}
\RequirePackage{grffile}
\RequirePackage[all]{xy}
\RequirePackage{wrapfig}
\RequirePackage{fancyvrb}
\RequirePackage{listings}

% Conditionals
\RequirePackage{ifthen}

% Header & Page Setup
\RequirePackage{fancyhdr}
\ifthenelse{\equal{\@oplargemargins}{1}}{}{\RequirePackage{fullpage}}

% Links
\RequirePackage{hyperref}

% ----- Questions -------------------------------------------------------------
\newcounter{questionCounter}
\newcounter{partCounter}[questionCounter]

% Prefix for questions
\newcommand{\questiontype}[0]{Question}

% Use this if your "written" questions are all under one section
% For example, if the homework handout has Section 5: Written Questions
% and all questions are 5.1, 5.2, 5.3, etc. set this to 5
% Use for 0 no prefix. Redefine as needed per-question.
\newcommand{\writtensection}[0]{0}

% Dispatch to the numbered or starred question command.
\newcommand{\question}{\@ifstar{\question@star}{\question@numbered}}

% Numbered question
\newcommand{\question@numbered}{%
  % Emit \newpage if option `newpage` is present
  \ifthenelse{\equal{\@opnewpage}{1}}{%
    \newpage
  }{}

  % Wrap in minipage so that we don't get a line break enywhere in between
  \begin{minipage}{\linewidth}%
    \stepcounter{questionCounter}%
      \vspace{.2in}%
      \ifx\writtensection\undefined{}
        \noindent{\bf \questiontype\ \arabic{questionCounter}. }%
        \else
          \ifnum\writtensection=0
          \noindent{\bf \questiontype\ \arabic{questionCounter}.}%
          \else
          \noindent{\bf \questiontype\ \writtensection.\arabic{questionCounter}}%
        \fi
      \vspace{0.3em} \hrule \vspace{.1in}%
  \end{minipage}
}

% Named question, takes one argument
\newcommand{\question@star}[1]{%
  % Emit \newpage if option `newpage` is present
  \ifthenelse{\equal{\@opnewpage}{1}}{%
    \newpage%
  }{}%
  % Wrap in minipage so that we don't get a line break enywhere in between
  \begin{minipage}{\linewidth}%
    \addtocounter{questionCounter}{1}%
    \setcounter{partCounter}{0}%
    \vspace{.2in}%
    \noindent{\bf \questiontype\ \arabic{questionCounter}.\hfill [#1 pts.]}%
    \vspace{0.3em} \hrule \vspace{.1in}%
  \end{minipage}
}

% Override normal section defintions
\renewcommand{\section}{\@ifstar{\section@star}{\question}}
\newcommand{\section@star}[1]{\question*{#1}}

% ----- Question Parts --------------------------------------------------------

\newenvironment{alphaparts}[0]{%
  \begin{enumerate}[label=\textbf{(\alph{partCounter})}]%
}{\end{enumerate}}

\newenvironment{arabicparts}[0]{%
  \begin{enumerate}[label=\textbf{\arabic{questionCounter}.\arabic{partCounter}})]%
}{\end{enumerate}}

\newcommand{\questionpart}[0]{\stepcounter{partCounter}\item}

% ----- Induction Environment -------------------------------------------------

\newenvironment{induction}[0]{%
  \begin{description}
}{\end{description}}

\newcommand{\basecase}{\item[Base Case]\mbox{}\\}
\newcommand{\indhyp}{\item[Induction Hypothesis]\mbox{}\\}
\newcommand{\indstep}{\item[Induction Step]\mbox{}\\}

% ----- Answer Box ------------------------------------------------------------

\newcommand{\answerbox}[1]{%
\begin{framed}
\vspace{#1}
\end{framed}}

% ----- Page Setup ------------------------------------------------------------

% Use block style paragraphs
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt plus 1pt}

\def\indented#1{\list{}{}\item[]}
\let\indented=\endlist

% ----- Title & Header --------------------------------------------------------
\pagestyle{empty}
\pagestyle{fancy}

%\if\@opanon%
\ifthenelse{\equal{\@opanon}{0}}{%
  \renewcommand{\maketitle}[0]{%
    % Setup header
    \setlength{\headheight}{15.2pt}
    \setlength{\headsep}{0.2in}
    \lhead{\hwclass{} }%
    \chead{\hwname{} }%
    \rhead{\hwtype{} \hwnum}%

    % Setup hrule in header
    \renewcommand{\headrulewidth}{0pt}
    \headrule{}

    % Don't put header on first page
    \thispagestyle{plain}

    \begin{center}
      {\Large \hwclass{} \hwtype{} \hwnum}

      \hwname{} 

      Due on \duedate
    \end{center}
    \renewcommand{\headrulewidth}{0.4pt}
  }

}%
{%
  \renewcommand{\maketitle}[0]{%
    % Make all pages plain
    \pagestyle{plain}

    % Put header on it's own page
    \begin{center}
      {\Large \hwclass{} \hwtype{} \hwnum}

      \hwname{} (\hwemail)

      Due on \duedate
    \end{center}
    \renewcommand{\headrulewidth}{0.4pt}
    \newpage
  }
}

% ----- For usage with pandoc converted documents -----------------------------

\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}

% -----------------------------------------------------------------------------
