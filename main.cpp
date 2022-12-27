# include <QtWidgets/QApplication>
# include <QtWidgets/QMainWindow>
# include <QtWidgets/QVBoxLayout>
//# include <QtWidgets/QAction>
# include <QtWidgets/QMenuBar>
# include <QtWidgets/QMessageBox>
# include "general.h"
# include "thread.h"

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);
  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  general_data *arg = new general_data;
  General *graph_area = new General (window);
  graph_area->set_arg (arg);
  QAction *action;

  if (graph_area->parse_command_line (argc, argv))
    {
      /*QMessageBox::warning (0, "Wrong input arguments!",
                            "Wrong input arguments!");*/
      return -1;
    }

  action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("0"));

  action = tool_bar->addAction ("&Change f->Pf", graph_area, SLOT (change_method ()));
  action->setShortcut (QString ("1"));

  action = tool_bar->addAction ("&Zoom_in", graph_area, SLOT (zoom ()));
  action->setShortcut (QString ("2"));

  action = tool_bar->addAction ("&Zoom_out", graph_area, SLOT (zoom_out ()));
  action->setShortcut (QString ("3"));

  action = tool_bar->addAction ("&n *= 2", graph_area, SLOT (increase_nx_and_ny ()));
  action->setShortcut (QString ("4"));

  action = tool_bar->addAction ("&n /= 2", graph_area, SLOT (reduce_nx_and_ny ()));
  action->setShortcut (QString ("5"));

  action = tool_bar->addAction ("&f + p*max", graph_area, SLOT (increase_splash ()));
  action->setShortcut (QString ("6"));

  action = tool_bar->addAction ("&f - p*max", graph_area, SLOT (reduce_splash ()));
  action->setShortcut (QString ("7"));

  action = tool_bar->addAction ("&left_rotate", graph_area, SLOT (left_rotate ()));
  action->setShortcut (QString ("8"));

  action = tool_bar->addAction ("&right_rotate", graph_area, SLOT (right_rotate ()));
  action->setShortcut (QString ("9"));

  action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->resize (800, 600);
  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  int p = graph_area->get_p ();
  pthread_t *tids = new pthread_t [p];
  args *arg_for_thread = new args [p];
  int nx = graph_area->get_nx ();
  int ny = graph_area->get_ny ();

  general_data *a = graph_area->get_general_data ();

  int i = 0;
  for (i = 0; i < p; i ++)
    {
      arg_for_thread [i].set_arg (p, i, nx, ny, graph_area->get_eps(), a);
    }
  for (int i = 0; i < graph_area->get_p (); i ++)
    {
      pthread_create(tids + i, 0, (&thread_func) , arg_for_thread + i);
    }

  QTimer *timer = new QTimer (graph_area);
  timer->setInterval (100);
  window->connect (timer, SIGNAL(timeout ()), graph_area, SLOT (timer ()));
  timer->start ();
  window->show ();

  app.exec ();
  delete [] tids;
  delete [] arg_for_thread;
  delete arg;
  delete window;
  return 0;
}
