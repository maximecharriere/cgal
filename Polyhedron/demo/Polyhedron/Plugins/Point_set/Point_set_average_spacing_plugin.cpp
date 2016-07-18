#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>
// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

using namespace CGAL::Three;
class Polyhedron_demo_point_set_average_spacing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionAverageSpacing;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionAverageSpacing = new QAction(tr("Average Spacing"), mainWindow);
    actionAverageSpacing->setProperty("subMenuName","Point Set Processing");
    actionAverageSpacing->setObjectName("actionAverageSpacing");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionAverageSpacing;
  }
  
  //! Applicable if the currently selected item is a
  //! points_with_normal_item.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionAverageSpacing_triggered();
}; // end Polyhedron_demo_point_set_average_spacing_plugin

void Polyhedron_demo_point_set_average_spacing_plugin::on_actionAverageSpacing_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    bool ok;

    const int nb_neighbors =
      QInputDialog::getInt((QWidget*)mw,
                               tr("Average Spacing"), // dialog title
                               tr("Number of neighbors:"), // field label
                               6, // default value = 1 ring
                               6, // min
                               1000, // max
                               1, // step
                               &ok);
    if(!ok) 
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Average spacing (k=" << nb_neighbors <<")...\n";

    Point_set::iterator points_begin = (points->nb_selected_points() == 0
                                        ? points->begin() : points->first_selected());

    // Computes average spacing
    double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
                                      points_begin, points->end(),
                                      nb_neighbors);

    // Print result
    Kernel::Sphere_3 bsphere = points->bounding_sphere();
    double radius = std::sqrt(bsphere.squared_radius());
    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Average spacing = " << average_spacing 
                                      << " = " << average_spacing/radius << " * point set radius (" 
                                      << task_timer.time() << " seconds, "
                                      << (memory>>20) << " Mb allocated)"
                                      << std::endl;
    QApplication::restoreOverrideCursor();

    QMessageBox::information(NULL,
                             tr("Average Spacing"),
                             tr("Average Spacing = %1 = %2 * point set radius")
                             .arg(average_spacing)
                             .arg(average_spacing/radius));
  }
}


#include "Point_set_average_spacing_plugin.moc"
